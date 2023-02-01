#version 440

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform float instrument_altitude;
uniform float instrument_area;

uniform float map_delta_az;
uniform bool  is_point_source;

uniform int   atm_nb_altitudes;
uniform int   atm_nb_angles;

uniform bool  use_aerosol;

uniform float src_altitude;

/* uniform vec3  instrument_los; */

/* uniform float instrument_fov; */
/* uniform float wavelength; */
/* uniform int   los_nb_points; */
/* uniform float instrument_lon; */
/* uniform float instrument_lat; */

#define X %%local_sizeX%%
#define Y %%local_sizeY%%
#define Z %%local_sizeZ%%

// Set up our compute groups
layout(local_size_x = X, local_size_y = Y, local_size_z = Z) in;

// Define constant values
const float EARTH_RADIUS = 6371.; //km
const float PI = 3.141592653589793;
const float DtoR = PI / 180;
const float RtoD = 1. / DtoR;


////////////////////////////////////////////////////////////////////////////////
// Structures of the input buffers. On for each array length.
////////////////////////////////////////////////////////////////////////////////

struct ScatteringData
{
    float angles;
    float aer_Pfct;
    float aer_Pfct_DoLP;
};

struct LOSData
{
    float altitudes;
    float volumes;
    float transmittance;
};

struct AtmosphereData
{
    float altitudes;
    float total_abs;
    float ray_beta;
    float aer_beta;
};

struct MapData
{
    float azimuts;
    float elevations;
};


////////////////////////////////////////////////////////////////////////////////
// Buffers for input and output data. Bindings must match the order of the buffer list of the wraper (all inputs first, then outputs).
////////////////////////////////////////////////////////////////////////////////

layout(std430, binding=0) buffer V_data_in{
    float data[];
} emission_data;

layout(std430, binding=1) buffer sca_data_in{
    ScatteringData data[];
} sca_data;

layout(std430, binding=2) buffer atm_data_in{
    AtmosphereData data[];
} atm_data;

layout(std430, binding=3) buffer los_data_in{
    LOSData data[];
} los_data;

layout(std430, binding=4) buffer map_elev_in{
    float data[];
} elev_data;

layout(std430, binding=5) buffer map_azim_in{
    float data[];
} az_data;

layout(std430, binding=6) buffer V_data_out{
    float data[];
} observation_data;


////////////////////////////////////////////////////////////////////////////////
// Function definitions
////////////////////////////////////////////////////////////////////////////////

float GetRange(float instrument_altitude, float instrument_elevation, float sca_altitude){
  float h = sca_altitude - instrument_altitude;
  return h / sin(instrument_elevation);
}

float GetAoLP(float src_a, float src_e, float i_a, float i_e){
  // Compute the AoLP of rayleigh scattering for a source at azimut src_a and elevation src_e as seen form the instrument. And an instrument pointing in the direction i_a, i_e
  vec3 src_uen = vec3(  sin(src_e),
                        cos(src_e) * sin(src_a),
                        cos(src_e) * cos(src_a));
  float Ce = cos(i_e);
  float Se = sin(i_e);
  float Ca = cos(i_a);
  float Sa = sin(i_a);
  float sin_AoRD = Se*(src_uen.z*Ca + src_uen.y*Sa) - src_uen.x*Ce;
  float cos_AoRD = src_uen.z*Sa - src_uen.y*Ca;

  return atan(sin_AoRD, cos_AoRD);
}

float GetSrcDistance(float src_elevation, float src_altitude, float instrument_altitude){
  float dist = - (EARTH_RADIUS + instrument_altitude) * sin(src_elevation);
  float square = - cos(src_elevation)*cos(src_elevation) * (EARTH_RADIUS + instrument_altitude)*(EARTH_RADIUS + instrument_altitude);
  square += (EARTH_RADIUS + src_altitude) * (EARTH_RADIUS + src_altitude);
  dist += sqrt(square);
  return  dist;
}


float GetEmissionArea(int ielev, float src_alt, float instrument_altitude){
  //Return the area of a pixel on the sky emission map in m**2.

  if(is_point_source){
    // return vec3(1, 1, 1);
    return 1;
  }

	float emin = elev_data.data[ielev];
  float emax = elev_data.data[ielev + 1];

  // Get the angle at the earth center between the instrument position and the source.
  // np.arcsin(o.AH_norm * np.cos(e) / (RT + h))
  float true_emin = asin(GetSrcDistance(emin, src_alt, instrument_altitude) * cos(emin) / (EARTH_RADIUS + src_altitude));
  float true_emax = asin(GetSrcDistance(emax, src_alt, instrument_altitude) * cos(emax) / (EARTH_RADIUS + src_altitude));

	float area = (EARTH_RADIUS + src_altitude)*(EARTH_RADIUS + src_altitude) * (cos(true_emax) - cos(true_emin)) * map_delta_az;

	return area * 1000000; //in m**2
}

int GetScaAngleIndex(float angle){
  for(int i = 0; i < atm_nb_angles - 1; i+=1){
    float mid_alt = (sca_data.data[i].angles + sca_data.data[i+1].angles) / 2;

    if (angle <= mid_alt){
      return i;
    }
  }
  return atm_nb_angles - 1;
}

int GetAltitudeIndex(float alt){
  //Returns the index of the altitude just below the input altitude
  for(int i = 0; i < atm_nb_altitudes - 1; i+=1){
    // float mid_alt = (atm_data.data[i].altitudes + atm_data.data[i+1].altitudes) / 2.;
    float mid_alt = atm_data.data[i].altitudes;

    if (alt <= mid_alt){
      return i-1;
    }
  }
  return atm_nb_altitudes - 1;
}

float GetInterpolationCoeff(float value, float lim1, float lim2){
  if(lim1 != lim2){
    return  (value - lim1) / (lim2 - lim1);
  }
  else{
    return 1.;
  }
}

float GetAtmosphereAbsorption(float alt1, float alt2){

  int ind1 = GetAltitudeIndex(alt1);
  int ind2 = GetAltitudeIndex(alt2);
  int ind11 = ind1+1;
  int ind21 = ind2+1;
  
  if(ind1 == atm_nb_altitudes - 1){ind11 -= 1;}
  if(ind2 == atm_nb_altitudes - 1){ind21 -= 1;}

  float abs1 = mix(atm_data.data[ind1].total_abs, atm_data.data[ind11].total_abs, GetInterpolationCoeff(alt1, atm_data.data[ind1].altitudes, atm_data.data[ind11].altitudes));
  float abs2 = mix(atm_data.data[ind2].total_abs, atm_data.data[ind21].total_abs, GetInterpolationCoeff(alt2, atm_data.data[ind2].altitudes, atm_data.data[ind21].altitudes));
  
  return abs(abs1 - abs2);
}

vec3 GetVParamFromLightParam(float I0, float DoLP, float AoLP){
  // Returns  the V parameters for a given radiant flux (any unit), a DoLP (between 0 and 1) and an AoLP in radians
  float V = I0 / 2.;
  float Vcos = I0 * DoLP * cos(2 * AoLP) / 4.;
  float Vsin = I0 * DoLP * sin(2 * AoLP) / 4.;
  return vec3(V, Vcos, Vsin);
}

float RSPhaseFunctionDoLP(float RD_angle){
  float Sa = sin(RD_angle);
  float Ca = cos(RD_angle);
  return Sa*Sa / (1 + Ca*Ca);
}

float RSPhaseFunction(float theta){
  // Simple approx
  float A = 3. / 4;
  float B = 1 + cos(theta)*cos(theta);

  // Chandrasekhar formula
  // if wavelength == 391.4:
  // 	gamma = 1.499 * 0.01
  // elif wavelength == 427.8:
  // 	gamma = 1.483 * 0.01
  // elif wavelength in [557.7, 550, 500]:
  // 	gamma = 1.442 * 0.01
  // elif wavelength == 630:
  // 	gamma = 1.413 * 0.01
  // A = 3 / (4 + 8 * gamma)
  // B = (1 + 3 * gamma) + (1 - gamma) * np.cos(theta) ** 2

  return A * B;
}

////////////////////////////////////////////////////////////////////////////////
// Main function
////////////////////////////////////////////////////////////////////////////////


void main()
{
/*
gl_NumWorkGroups
    This variable contains the number of work groups passed to the dispatch function.
gl_WorkGroupID
    This is the current work group for this shader invocation. Each of the XYZ components will be on the half-open range [0, gl_NumWorkGroups.XYZ).
gl_LocalInvocationID
    This is the current invocation of the shader within the work group. Each of the XYZ components will be on the half-open range [0, gl_WorkGroupSize.XYZ).
gl_GlobalInvocationID
    This value uniquely identifies this particular invocation of the compute shader among all invocations of this compute dispatch call. It's a short-hand for the math computation:
        gl_WorkGroupID * gl_WorkGroupSize + gl_LocalInvocationID;
gl_LocalInvocationIndex
    This is a 1D version of gl_LocalInvocationID. It identifies this invocation's index within the work group. It is short-hand for this math computation:

      gl_LocalInvocationIndex =
          gl_LocalInvocationID.z * gl_WorkGroupSize.x * gl_WorkGroupSize.y +
          gl_LocalInvocationID.y * gl_WorkGroupSize.x +
          gl_LocalInvocationID.x;

Shader wrap creates a 2d grid of WorkGroups of size (N_elevations, N_azimuts). Each work group contains as many invocation as there are points along the line of sight of the instrument.
 */
    
    


    // Setting and computing indices and correspoinding paramters
    int src_elevation_index = int(gl_GlobalInvocationID.x);
    float src_elevation = (elev_data.data[src_elevation_index+1] + elev_data.data[src_elevation_index])/2.;

    int src_azimut_index = int(gl_GlobalInvocationID.y);
    float src_azimut = (az_data.data[src_azimut_index+1] + az_data.data[src_azimut_index])/2.;


   float src_distance = GetSrcDistance(src_elevation, src_altitude, instrument_altitude);
   float src_horiz_distance = sqrt(src_distance*src_distance - src_altitude*src_altitude);
  //  float src_horiz_distance = src_altitude / tan(src_elevation);
  //  float src_distance = src_altitude / sin(src_elevation);

    float azimut_diff = src_azimut - instrument_azimut;
    int sca_index = int(gl_LocalInvocationID.z);
    float sca_altitude = los_data.data[sca_index].altitudes;
    if(sca_altitude == 0){return;}
    int sca_altitude_index = GetAltitudeIndex(sca_altitude);
    float sca_volume = los_data.data[sca_index].volumes;
    float sca_range = GetRange(instrument_altitude, instrument_elevation, sca_altitude);
    if(sca_range == 0){return;}
    float sca_range_horiz = sca_range * cos(instrument_elevation);


    int out_buffer_index = 3 * int(src_elevation_index * gl_NumWorkGroups.y * gl_WorkGroupSize.z + src_azimut_index * gl_WorkGroupSize.z + sca_index);
    int V_index    = out_buffer_index;
    int Vcos_index = out_buffer_index + 1;
    int Vsin_index = out_buffer_index + 2;

    int in_buffer_index = 3 * int(src_elevation_index * gl_NumWorkGroups.y + src_azimut_index);

    // From now on, the instrument is at point A, the scatering volume is centered on point S, the emission on point E. So elevations and angles use these points name.

// Doing some geometry computations
    float SE_horiz_square = src_horiz_distance*src_horiz_distance + sca_range_horiz*sca_range_horiz - 2*src_horiz_distance*sca_range_horiz*cos(azimut_diff);
    //SE is doesn't exactly match the CPU value for approx reasons in the geometry made here. Could be corrected.
    float SE = sqrt(SE_horiz_square + (src_altitude - sca_altitude)*(src_altitude - sca_altitude));
    if(SE == 0){return;}

    float sca_angle = (sca_range * sca_range + SE * SE - src_distance * src_distance) / (2 * SE * sca_range);
    sca_angle = PI - acos(sca_angle);
    int sca_angle_index = GetScaAngleIndex(sca_angle);

// Get AoLP from previous geometry.
// WARNING: the AoLP is not exactly identical than with CPU function (<1deg difference). To solve, check geometry simplifications (flat earth,etc...)
    float AoLP = GetAoLP(src_azimut, src_elevation, instrument_azimut, instrument_elevation);

    ////////////// GOOD //////////////////


// Propagate light from emission to instrument
    // init the source radiance
    float I0 = emission_data.data[in_buffer_index]; //[nW/m2/sr]

    // Normalize by emission area
    I0 *= GetEmissionArea(src_elevation_index, src_altitude, instrument_altitude); // [nW / m2 / sr * m2] = [nW / sr]

    // elevation square law to scattering point
    if (SE > 0){
      I0 /= (SE * SE); // [nW / m2 / sr * m2 / km2] = [nW / km2]
    }
		else{
			I0 = 0;
    }

    // float debug2 = GetEmissionArea(src_elevation_index, src_altitude, instrument_altitude).z;

    // Losses between emission and scattering
    float opt_depth = 0;
    if (sca_altitude != 0){
      float delta_z = abs(src_altitude - sca_altitude);
      opt_depth = GetAtmosphereAbsorption(src_altitude, sca_altitude) * SE / delta_z;
    }

    // float debug1 = -opt_depth;

    I0 *= exp(-opt_depth); // [nW / km2]

    if (sca_range != 0){
      I0 *= sca_volume * instrument_area / (sca_range*sca_range) / 4 / PI;
    }

    // Scattering parameters (phase function and cross section) for Rayleigh and aerosols
		float	Crs = atm_data.data[sca_altitude_index].ray_beta; //in km-1
    /* float	ray_Pfct = atm_ray_Pfct[sca_angle_index]; // in sr */
    float	ray_Pfct = RSPhaseFunction(sca_angle); // in sr
    /* float	ray_Pfct_DoLP = atm_ray_Pfct_DoLP[sca_angle_index]; */
    float	ray_Pfct_DoLP = RSPhaseFunctionDoLP(sca_angle);


		float Caer = 0;
    float Paer = 0;
    float DoLP_aer = 0;
		if (use_aerosol){
			Caer     = atm_data.data[sca_altitude_index].aer_beta; //in km-1
			Paer     = sca_data.data[sca_angle_index].aer_Pfct;
			DoLP_aer = sca_data.data[sca_angle_index].aer_Pfct_DoLP;
    }

		float I0_rs  = I0 * Crs  * ray_Pfct;  // [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]
		float I0_aer = I0 * Caer * Paer;      // [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]

    I0 = I0_rs + I0_aer;

    I0 *= los_data.data[sca_index].transmittance;


  	// Mix of flux and DoLP when using aerosols
    float DoLP = 0;
    if (I0 != 0){
      DoLP = (I0_rs * ray_Pfct_DoLP + I0_aer * DoLP_aer) / (I0_rs + I0_aer);
    }

    // float debug1 = DoLP;
    // For Mie scattering, a DoLP < 0 = AoLP is parallel to scatering plane ! So 90 degrees from Rayleigh scattering.
    if (DoLP < 0){
      AoLP += PI/2;
      DoLP *= -1;
    }



		vec3 stokes_param = GetVParamFromLightParam(I0, DoLP, AoLP);

    observation_data.data[V_index]    = stokes_param.x;
    observation_data.data[Vcos_index] = stokes_param.y;
    observation_data.data[Vsin_index] = stokes_param.z;

    // observation_data.data[V_index]    = debug1;
    // observation_data.data[Vcos_index] = debug2;
}

























































