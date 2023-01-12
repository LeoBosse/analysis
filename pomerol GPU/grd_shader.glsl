#version 440

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform float instrument_altitude;
uniform int   atm_nb_altitudes;
uniform int   atm_nb_angles;

uniform float instrument_area;
uniform float map_delta_az;
uniform bool  is_point_source;
uniform bool  use_aerosol;
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
const float PI = 3.14159;
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
    float distances;
};


////////////////////////////////////////////////////////////////////////////////
// Buffers for input and output data. Bindings must match the order of the buffer list of the wraper.
////////////////////////////////////////////////////////////////////////////////

layout(std430, binding=0) buffer V_data_in{
    float data[];
} emission_data;

layout(std430, binding=1) buffer V_data_out{
    float data[];
} observation_data;

layout(std430, binding=2) buffer sca_data_in{
    ScatteringData data[];
} sca_data;

layout(std430, binding=3) buffer atm_data_in{
    AtmosphereData data[];
} atm_data;

layout(std430, binding=4) buffer los_data_in{
    LOSData data[];
} los_data;

layout(std430, binding=5) buffer map_dist_in{
    float data[];
} dist_data;

layout(std430, binding=6) buffer map_azim_in{
    float data[];
} az_data;

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

float GetEmissionArea(int idist){
  //Return the area of a pixel on the ground emission map in m**2.
  if (is_point_source){
    return 1;
  }

  float area = 0.5 * abs(map_delta_az);
  area *= (dist_data.data[idist+1]*dist_data.data[idist+1] - dist_data.data[idist]*dist_data.data[idist]);
  area *= 1000000;

  return area;
}

int GetScaAngleIndex(float angle){
  for(int i = 0; i < atm_nb_angles - 1; i++){
    if (angle <= (sca_data.data[i].angles + sca_data.data[i+1].angles) / 2){
      return i;
    }
  }
  return atm_nb_angles;
}

int GetAltitudeIndex(float alt){
  for(int i = 0; i < atm_nb_altitudes - 1; i+=1){
    float mid_alt = (atm_data.data[i].altitudes + atm_data.data[i+1].altitudes) / 2.;

    if (alt <= mid_alt){
      return i;
    }
  }
  return atm_nb_altitudes;
}

float GetAtmosphereAbsorption(float alt1, float alt2){

  float abs1 = atm_data.data[GetAltitudeIndex(alt1)].total_abs;
  float abs2 = atm_data.data[GetAltitudeIndex(alt2)].total_abs;

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

Shader wrap creates a 2d grid of WorkGroups of size (N_distances, N_azimuts). Each work group contains as many invocation as there are points along the line of sight of the instrument.
 */
    // Setting and computing indices and correspoinding paramters
    int src_distance_index = int(gl_GlobalInvocationID.x);
    float src_distance = (dist_data.data[src_distance_index+1] + dist_data.data[src_distance_index])/2.;

    int src_azimut_index = int(gl_GlobalInvocationID.y);
    float src_azimut = (az_data.data[src_azimut_index+1] + az_data.data[src_azimut_index])/2.;
    float src_altitude = 0.;

    float azimut_diff = src_azimut - instrument_azimut;
    int sca_index = int(gl_LocalInvocationID.z);
    float sca_altitude = los_data.data[sca_index].altitudes;
    if(sca_altitude == 0){return;}
    int sca_altitude_index = GetAltitudeIndex(sca_altitude);
    float sca_volume = los_data.data[sca_index].volumes;
    float sca_range = GetRange(instrument_altitude, instrument_elevation, sca_altitude);
    if(sca_range == 0){return;}
    float sca_range_horiz = sca_range * cos(instrument_elevation);


    int out_buffer_index = 3 * int(src_distance_index * gl_NumWorkGroups.y * gl_WorkGroupSize.z + src_azimut_index * gl_WorkGroupSize.z + sca_index);
    int V_index    = out_buffer_index;
    int Vcos_index = out_buffer_index + 1;
    int Vsin_index = out_buffer_index + 2;

    int in_buffer_index = 3 * int(src_distance_index * gl_NumWorkGroups.y + src_azimut_index);

    // From now on, the instrument is at point A, the scatering volume is centered on point S, the emission on point E. So distances and angles use these points name.

// Doing some geometry computations
    float SE_horiz_square = src_distance*src_distance + sca_range_horiz*sca_range_horiz - 2*src_distance*sca_range_horiz*cos(azimut_diff);
    //SE is doesn't exactly match the CPU value for approx reasons in the geometry made here. Could be corrected.
    float SE = sqrt(SE_horiz_square + sca_altitude*sca_altitude);
    if(SE == 0){return;}

    float sca_angle = (sca_range * sca_range + SE * SE - src_distance * src_distance) / (2 * SE * sca_range);
    sca_angle =  acos(sca_angle);

    int sca_angle_index = GetScaAngleIndex(sca_angle);
    float SAE = (sca_range * sca_range + src_distance * src_distance - SE * SE) / (2 * src_distance * sca_range);
    SAE = acos(SAE);

    /* vec3 sca_plane_normal = vec3( instrument_los[1] * sin(azimut_diff) - instrument_los[2] * cos(azimut_diff),
    -instrument_los[0] * sin(azimut_diff),
    instrument_los[0] * cos(azimut_diff)); */



// Get AoLP from previous geometry.
// WARNING: the AoLP is not exactly identical than with CPU function (<1deg difference). To solve, check geometry simplifications (flat earth,etc...)
    float AoLP = GetAoLP(src_azimut, 0., instrument_azimut, instrument_elevation);

// Propagate light from emission to instrument
    // init the source radiance
    float I0 = emission_data.data[in_buffer_index]; //[nW/m2/sr]


    // Normalize by emission area
    I0 *= GetEmissionArea(src_distance_index); // [nW / m2 / sr * m2] = [nW / sr]


    // distance square law to scattering point
    if (SE > 0){
      I0 /= (SE * SE); // [nW / m2 / sr * m2 / km2] = [nW / km2]
    }
		else{
			I0 = 0;
    }

    // Losses between emission and scattering
    float opt_depth = 0;
    if (sca_altitude != 0){
      float delta_z = abs(src_altitude - sca_altitude);
      opt_depth = GetAtmosphereAbsorption(src_altitude, sca_altitude) * SE / delta_z;
    }

    I0 *= exp(- opt_depth); // [nW / km2]

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
			Caer = atm_data.data[sca_altitude_index].aer_beta; //in km-1
			Paer = sca_data.data[sca_angle_index].aer_Pfct;
			DoLP_aer = sca_data.data[sca_angle_index].aer_Pfct_DoLP;
    }

		float I0_rs  = I0 * Crs  * ray_Pfct;  // [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]
		float I0_aer = I0 * Caer * Paer;      // [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]

    I0 = I0_rs + I0_aer;

    I0 *= los_data.data[sca_index].transmittance;
    ////////////// GOOD //////////////////

  	// Mix of flux and DoLP when using aerosols
    float DoLP = 0;
    if (I0 != 0){
      DoLP = (I0_rs * ray_Pfct_DoLP + I0_aer * DoLP_aer) / (I0_rs + I0_aer);
    }

    // For Mie scattering, a DoLP < 0 = AoLP is parallel to scatering plane ! So 90 degrees from Rayleigh scattering.
    if (DoLP < 0){
      AoLP += PI/2;
      DoLP *= -1;
    }

    float debug1 = DoLP;

		vec3 stokes_param = GetVParamFromLightParam(I0, DoLP, AoLP);

    observation_data.data[V_index]    = stokes_param.x;
    observation_data.data[Vcos_index] = stokes_param.y;
    observation_data.data[Vsin_index] = stokes_param.z;

    /* observation_data.data[V_index]    = debug1; */
    /* observation_data.data[Vcos_index] = debug2; */
}




////////////////////////////////////////////////////////////////////////////////
// Main function as done previously on the CPU
////////////////////////////////////////////////////////////////////////////////

/*
		I0 = self.ground_map.cube[time, idist, iaz] #[nW/m2/sr]
		alt = self.altitudes[ialt]
  	AR, RE, RD_angle,RAE, alt_E = self.GetGeometryFromAzDist(a_E, e_E, alt, ret_RAE=True)
    vol = self.atmosphere.GetVolume(AR, self.ouv_pc, index = ialt, unit="km", type="index") #km3
    AoLP = obs.GetRayleighAngle(a_E, e_E)

    obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, self.a_pc, self.e_pc, A_alt = self.instrument_altitude, init_full=False)

    alt = obs.GetPCoordinatesFromRange(AR)[2]
    I0 *= self.ground_map.GetArea(idist) # [nW / m2 / sr * m2] = [nW / sr]

    if RE > 0:
    I0 /= RE ** 2 # [nW / m2 / sr * m2 / km2] = [nW / km2]
    else:
    I0 = 0

    opt_depth = 0
    O3_abs = 0
    aer_abs = 0
    if alt != 0:
    delta_z = abs(alt_E - alt)
    opt_depth = self.atmosphere.GetRSOpticalDepth(alt_E, alt) * RE / delta_z
    if self.atmosphere.use_ozone:
    O3_abs = self.atmosphere.GetO3Absorbtion(alt_E, alt) * RE / delta_z
    if self.atmosphere.use_aerosol:
    aer_abs = self.atmosphere.GetAerosolsAbsorbtion(alt_E, alt) * RE / delta_z
    # print("DEBUG opt_depth: ER", opt_depth, 1-opt_depth, np.exp(-opt_depth))

    I0 *= np.exp(- opt_depth - O3_abs - aer_abs) # [nW / km2]


    Crs = self.atmosphere.GetRSVolumeCS(alt) #in km-1
    P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)# in sr

    Caer, Paer, DoLP_aer = 0, 0, 0
    if self.atmosphere.use_aerosol:
    Caer = self.atmosphere.GetAerosolCS(alt) #in km-1
    Paer, DoLP_aer = self.atmosphere.GetAerosolPhaseFunction(RD_angle) # Paer in sr. pi=FRONTscattering, 0=backscattring

    I0_rs = 0
    I0_aer = 0
    if AR != 0:
    I0_aer = I0 * Caer * Paer * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi # [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]
    I0_rs = I0 * Crs * P * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi # [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]

    DoLP_rs = self.atmosphere.GetRSPhaseFunctionDoLP(RD_angle)

    I0 = I0_rs + I0_aer
    I0 *= self.atmosphere.los_transmittance[ialt]


    # Option 1: Mix of flux and DoLP when using aerosols
    if I0 != 0:
    DoLP = (I0_rs * DoLP_rs + I0_aer * DoLP_aer) / I0
    else:
    DoLP = 0

    if DoLP < 0: # For Mie scattering, a DoLP < 0 = AoLP is parallel to scatering plane ! So 90 degrees from Rayleigh scattering.
    AoLP += np.pi/2
    DoLP *= -1
    geo_list = [[self.a_pc, self.e_pc, AR, AR, RE, RD_angle, alt, AoLP, vol, self.atmosphere.d_los_list[ialt]]]
    for a_A, e_A, r, AR, RE, RD_angle, alt, AoLP, dvol, dr in geo_list:

    if RE == 0:
    print("WARNING: ER = 0 !!!")
    continue





    V, Vcos, Vsin = self.GetVParamFromLightParam(I0, DoLP, AoLP) #Using aerosols
    return V_T, Vcos_T, Vsin_T */
