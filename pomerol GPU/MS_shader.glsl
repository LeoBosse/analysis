#version 440

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform float instrument_altitude;
// uniform float instrument_area;

uniform int   atm_nb_altitudes;
// uniform int   atm_nb_angles;

uniform float min_altitude;
uniform float max_altitude;
// uniform float distance_limit;
uniform int   scattering_limit;

uniform float increment_length;

uniform int total_ray_number;

#define X %%local_sizeX%%
#define Y %%local_sizeY%%
#define Z %%local_sizeZ%%

#define HISTORY_MAX_SIZE %%scattering_limit%%

// Set up our compute groups
layout(local_size_x = X, local_size_y = Y, local_size_z = Z) in;

// Define constant values
const float EARTH_RADIUS = 6371.; //km
const float PI = 3.141592653589793;
const float DtoR = PI / 180;
const float RtoD = 1. / DtoR;

const int ESCAPE_IN_SPACE   = 0;
const int TOUCH_GROUND      = 1;
const int VALID             = 2;
const int INVALID           = 3;


vec3 initial_vec = vec3(sin(instrument_elevation), 
                        cos(instrument_elevation) * sin(instrument_azimut), 
                        cos(instrument_elevation) * cos(instrument_azimut));

uint ray_index = gl_LocalInvocationIndex;


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
} history[HISTORY_MAX_SIZE];


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

// layout(std430, binding=3) buffer history_data_out{
//     RayHistoryData data[];
// } history;



////////////////////////////////////////////////////////////////////////////////
// Function definitions
////////////////////////////////////////////////////////////////////////////////

uint _pcg(uint seed){
  uint state = seed * 747796405u + 2891336453u;
  uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

float RandomNumber(uint p){
	return float(_pcg(int(p))) / float(uint(0xffffffff));
}

float RandomNumber(vec2 p){
	return float(_pcg(_pcg(uint(p.x)) + uint(p.y))) / float(uint(0xffffffff));
}

vec3 MatrixProduct(mat3 M, vec3 v){
  M[0] *= v.x;
  M[1] *= v.y;
  M[2] *= v.z;
  return vec3(M[0].x+M[1].x+M[2].x , M[0].y+M[1].y+M[2].y , M[0].z+M[1].z+M[2].z);
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

// int GetScaAngleIndex(float angle){
//   for(int i = 0; i < atm_nb_angles - 1; i+=1){
//     float mid_alt = (sca_data.data[i].angles + sca_data.data[i+1].angles) / 2;

//     if (angle <= mid_alt){
//       return i;
//     }
//   }
//   return atm_nb_angles - 1;
// }

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

vec2 GetAtmosphereCrossSection(float alt){

  int ind = GetAltitudeIndex(alt);
  int ind1 = ind + 1;
  
  if(ind1 == atm_nb_altitudes - 1){ind1 -= 1;}

  float ray_cs = mix(atm_data.data[ind].ray_beta, atm_data.data[ind1].ray_beta, GetInterpolationCoeff(alt, atm_data.data[ind].altitudes, atm_data.data[ind1].altitudes));
  float aer_cs = mix(atm_data.data[ind].aer_beta, atm_data.data[ind1].aer_beta, GetInterpolationCoeff(alt, atm_data.data[ind].altitudes, atm_data.data[ind1].altitudes));
  
  return vec2(ray_cs, aer_cs);
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

vec4 GetScatteringPlane(){
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
    float plane_angle =  2 * PI * (0.5 - ray_index / total_ray_number);
    // Normal vector to plane in ref frame of instrument (==AoLP)
    vec3 normal = vec3(0, sin(plane_angle), cos(plane_angle));
    // float Pxi = 0;
    // float Pyi = sin(plane_angle);
    // float Pzi = cos(plane_angle);
    // Normal vector to plane in UpEastNorth coord.
    vec3 P_normal = MatrixProduct(Ria, normal);

    return vec4(P_normal, plane_angle);
}

void WriteHistory(int index, float altitude, vec3 vec, float sca_angle, float cs_ray, float cs_aer, float segment_length){
  history[index].altitude       = altitude;
  history[index].vec            = vec;
  history[index].sca_angle      = sca_angle;
  history[index].cs_ray         = cs_ray;
  history[index].cs_aer         = cs_aer;
  history[index].segment_length = segment_length;
}

void InitializeHistory(){
    for(int sca_ID=0; sca_ID < HISTORY_MAX_SIZE; sca_ID+=1){
      WriteHistory(sca_ID, instrument_altitude, initial_vec, -1, -1, -1, 0);
    }
}

bool IsScattered(float cross_section, int sca_ID){
   return RandomNumber(vec2(ray_index, sca_ID)) < cross_section;

}

float GetScatteringAngle(float cs_ray, float cs_aer){
    // float phase_function = cs_ray * self.atmosphere.profiles["ray_Phase_Fct"] + cs_aer * self.atmosphere.profiles["aer_Phase_Fct"];

    // phase_function /= (cs_ray + cs_aer)

    // proba = phase_function / sum(phase_function);

    // random_angle = np.random.choice([-1, 1]) * np.random.choice(self.atmosphere.profiles["sca_angle"], p = proba)
    return 0;
}

int Propagate(vec3 scattering_plane){
    
    bool is_scattered = false;
    vec2 cross_sections;

    int ray_state = VALID;

    int sca_ID = 0;
    float altitude  = history[0].altitude;
    vec3 vec        = history[0].vec;
    float segment_length = history[0].segment_length;

    // while the ray is valid == dont touch the ground or escape in the sky.
    while(ray_state == VALID && sca_ID < scattering_limit){
        // While the ray is not scattered == while the ray goes in a straight line == while we are on the same path segment
        while (ray_state == VALID && !is_scattered){

            altitude += (vec.x * increment_length);  //Increment the altitude
            segment_length += increment_length;      //Increment the segment length

            if (altitude < min_altitude){
                history[sca_ID] = RayHistoryData(altitude, vec, -1, -1, -1, segment_length);
                ray_state = TOUCH_GROUND;
            }
            if(altitude > max_altitude){
                ray_state = ESCAPE_IN_SPACE;
            }

            cross_sections = GetAtmosphereCrossSection(altitude) * increment_length;
            float total_cross_section = cross_sections.x + cross_sections.y;

            is_scattered = IsScattered(total_cross_section, sca_ID);
        }

        sca_ID += 1;
        float sca_angle = GetScatteringAngle(cross_sections.x, cross_sections.y);
        // vec3 vec = self.RotateAboutPlane(vec, scattering_plane, sca_angle)

        WriteHistory(sca_ID, altitude, vec, sca_angle, cross_sections.x, cross_sections.y, segment_length);
        
        segment_length = 0;

    }
    return ray_state;

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

    
    uint V_index    = ray_index;
    uint Vcos_index = ray_index + 1;
    uint Vsin_index = ray_index + 2;

    vec4 scattering_plane = GetScatteringPlane();
    vec3 scattering_plane_normal = scattering_plane.xyz;
    float AoLP = scattering_plane.w;

    bool invalid_ray = true;
    while (invalid_ray){
        InitializeHistory();
        Propagate(scattering_plane_normal);
    }

    vec3 stokes_param = vec3(1, 0, 0);

    observation_data.data[V_index]    = stokes_param.x;
    observation_data.data[Vcos_index] = stokes_param.y;
    observation_data.data[Vsin_index] = stokes_param.z;

    // observation_data.data[V_index]    = debug1;
    // observation_data.data[Vcos_index] = debug2;
}

























































