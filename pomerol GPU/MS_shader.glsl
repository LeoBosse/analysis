#version 440

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform float instrument_altitude;
uniform float instrument_area;

uniform float min_altitude;
uniform float max_altitude;
uniform float distance_limit;
uniform int   scattering_limit;

uniform float increment_length;

uniform int total_ray_number;

vec3 initial_vec = vec3(sin(instrument_elevation), 
                        cos(instrument_elevation) * sin(instrument_azimut), 
                        cos(instrument_elevation) * cos(instrument_azimut));


// uniform vec3 instrument_los;

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

const int ESCAPE_IN_SPACE   = 0;
const int TOUCH_GROUND      = 1;



////////////////////////////////////////////////////////////////////////////////
// Structures of the input buffers. On for each array length.
////////////////////////////////////////////////////////////////////////////////

struct ScatteringData
{
    float angles;
    float aer_Pfct;
    float aer_Pfct_DoLP;
};

struct AtmosphereData
{
    float altitudes;
    float total_abs;
    float ray_beta;
    float aer_beta;
};

struct RayHistoryData{
    float altitude;
    vec3  vec;
    float sca_angle;
    float cs_ray;
    float cs_aer;
    float segment_length;
} history[scattering_limit];


////////////////////////////////////////////////////////////////////////////////
// Buffers for input and output data. Bindings must match the order of the buffer list of the wraper (all inputs first, then outputs).
////////////////////////////////////////////////////////////////////////////////

layout(std430, binding=0) buffer sca_data_in{
    ScatteringData data[];
} sca_data;

layout(std430, binding=1) buffer atm_data_in{
    AtmosphereData data[];
} atm_data;


layout(std430, binding=2) buffer V_data_out{
    float data[];
} observation_data;


////////////////////////////////////////////////////////////////////////////////
// Function definitions
////////////////////////////////////////////////////////////////////////////////


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

vec4 GetScatteringPlane(int rayID){
    // Rotation matrix from instrument ref frame to UpEastNorth
    float ca = cos(instrument_azimut);
    float sa = sin(instrument_azimut);
    float ce = cos(instrument_elevation);
    float se = sin(instrument_elevation);

    mat3 Ria = mat3(se,     ce*sa,      ce*ca,      //First col
                    0,      -ca,        sa,         //Second col
                    ce,     -se*sa,     -se*ca);     //Third col
        // [	    [se, 		0,		ce],
        //          [ce * sa,	-ca,	-se * sa],
        //          [ce * ca,	sa,		-se * ca]])


    // Plane angle
    float plane_angle =  2 * PI * (0.5 - rayID / total_ray_number);
    // Normal vector to plane in ref frame of instrument (==AoLP)
    float Pxi = 0;
    float Pyi = sin(plane_angle);
    float Pzi = cos(plane_angle);
    // Normal vector to plane in UpEastNorth coord.
    float P_normal = np.dot(Ria, (Pxi, Pyi, Pzi));

    return vec4(P_normal, plane_angle);
}

void InitializeHistory(){
    for(int sca_ID; sca_ID < scattering_limit; sca_ID+=1){
        history[sca_ID].altitude = instrument_altitude;
        history[sca_ID].vec = initial_vec;
        history[sca_ID].sca_angle = -1;
        history[sca_ID].cs_ray = -1;
        history[sca_ID].cs_aer = -1;
        history[sca_ID].segment_length = -1;
    }
}

bool IsScattered(altitude){
    return 1;
}

int Propagate(int sca_ID, float altitude, vec3 vec, vec3 scattering_plane, float segment_length){
    
    bool is_scattered = 0;

    // While the ray is not scattered == while the ray goes in a straight line == while we are on the same path segment
    while (!is_scattered){

        altitude += (vec.x * increment_length)  //Increment the altitude
        segment_length += increment_length      //Increment the segment length

        if (altitude < min_altitude){
            hist[sca_ID] = RayHistoryData(altitude, vec, -1, -1, -1, segment_length);
            return TOUCH_GROUND
        }
        if(altitude > max_altitude){
            return ESCAPE_IN_SPACE
        }

        is_scattered = IsScattered(alt)
    }

    _, cs_ray, cs_aer = is_scattered
    sca_angle = self.GetScatteringAngle(cs_ray, cs_aer)
    # print(vec, sca_angle*RtoD)
    vec = self.RotateAboutPlane(vec, scattering_plane, sca_angle)
    # print(vec)

    hist = np.append(hist, [[alt, vec, sca_angle, cs_ray, cs_aer, segment_length]], axis = 0)

    # print(hist)

    segment_length = 0

    return self.Propagate(alt, vec, scattering_plane, id_ray=id_ray, hist = hist, segment_length=segment_length)



    return 0
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

    int ray_index = gl_LocalInvocationIndex;
    int V_index    = ray_index;
    int Vcos_index = ray_index + 1;
    int Vsin_index = ray_index + 2;


    vec4 scattering_plane = GetScatteringPlane(ray_index);
    vec3 scattering_plane_normal = scattering_plane.xyz;
    float AoLP = scattering_plane.w;

    bool invalid_ray = 1;
    while (invalid_ray){

        InitializeHistory();
        history = Propagate(scattering_plane_normal);


    }


    vec3 stokes_param = vec3(1, 0, 0);

    observation_data.data[V_index]    = stokes_param.x;
    observation_data.data[Vcos_index] = stokes_param.y;
    observation_data.data[Vsin_index] = stokes_param.z;

    // observation_data.data[V_index]    = debug1;
    // observation_data.data[Vcos_index] = debug2;
}

























































