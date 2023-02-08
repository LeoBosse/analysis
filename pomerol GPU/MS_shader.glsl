#version 440

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform float instrument_altitude;
// uniform float instrument_area;

uniform int   atm_nb_altitudes;
uniform float wavelength;
// uniform int   atm_nb_angles;

uniform float min_altitude;
uniform float max_altitude;
uniform float distance_limit;
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

const int INVALID           = -1;
const int VALID             = 1;
const int ESCAPE_IN_SPACE   = 100;
const int TOUCH_GROUND      = 0;


vec3 initial_vec = vec3(sin(instrument_elevation), 
                        cos(instrument_elevation) * sin(instrument_azimut), 
                        cos(instrument_elevation) * cos(instrument_azimut));

int nb_scattering_events = 0;

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
uint rays_per_work_group = gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z;
uint GlobalInvocationIndex = gl_WorkGroupID.z * gl_NumWorkGroups.x * gl_NumWorkGroups.y +
                             gl_WorkGroupID.y * gl_NumWorkGroups.x +
                             gl_WorkGroupID.x;
uint ray_index =  GlobalInvocationIndex * rays_per_work_group + gl_LocalInvocationIndex;


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

struct MapData
{
    float azimuts;
    float distances;
};

struct RayHistoryData{
    float altitude;
    vec3  vec;
    float sca_angle;
    float cs_ray;
    float cs_aer;
    float segment_length;
    vec3 sca_plane;
    int state;
    float plane_rotation_angle;
} history[HISTORY_MAX_SIZE];



////////////////////////////////////////////////////////////////////////////////
// Buffers for input and output data. Bindings must match the order of the buffer list of the wraper (all inputs first, then outputs, in orders).
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

layout(std430, binding=3) buffer debug_data_out{
    float data[];
} debug_data;


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

vec3 GetScatteringPlane(vec3 ray_vec, float plane_rotation_angle){
    // From an incident vector and an angle, returns the normal to the defined plane. The angle is the angle of the plane with the vertical (similar to the AoLP).

    // float ray_azimut = atan(ray_vec.y, ray_vec.z);
    float ray_elevation = asin(ray_vec.x);

    // Rotation matrix from ray ref frame (forward, horiz, vert) to UpEastNorth
    float ca = ray_vec.z; //cos(ray_azimut);
    float sa = ray_vec.y; //sin(ray_azimut);
    float ce = cos(ray_elevation);
    float se = ray_vec.x; //sin(ray_elevation);

    mat3 Ria = mat3(se,     ce*sa,      ce*ca,      //First col
                    0,      -ca,        sa,         //Second col
                    ce,     -se*sa,     -se*ca);    //Third col
        // [	    [se, 		0,		ce],
        //          [ce * sa,	-ca,	-se * sa],
        //          [ce * ca,	sa,		-se * ca]])


    // Normal vector to plane in ref frame of instrument (==AoLP)
    vec3 normal = vec3(0, sin(plane_rotation_angle), cos(plane_rotation_angle));
    // float Pxi = 0;
    // float Pyi = sin(plane_angle);
    // float Pzi = cos(plane_angle);
    // Normal vector to plane in UpEastNorth coord.
    vec3 P_normal = Ria * normal;

    return P_normal;
}

void WriteHistory(int index, float altitude, vec3 vec, float sca_angle, float cs_ray, float cs_aer, float segment_length, vec3 sca_plane, int ray_state, float plane_rotation_angle){
  history[index].altitude       = altitude;
  history[index].vec            = vec;
  history[index].sca_angle      = sca_angle;
  history[index].plane_rotation_angle = plane_rotation_angle;
  history[index].cs_ray         = cs_ray;
  history[index].cs_aer         = cs_aer;
  history[index].segment_length = segment_length;
  history[index].state          = ray_state;
}

void InitializeHistory(){
    for(int sca_ID=0; sca_ID < HISTORY_MAX_SIZE; sca_ID+=1){
      WriteHistory(sca_ID, instrument_altitude, initial_vec, -1, -1, -1, 0, vec3(0, 0, 0), VALID, 0);
    }
}

bool IsScattered(float cross_section, int sca_ID){
   return RandomNumber(vec2(ray_index, sca_ID)) < cross_section;

}

vec3 RotateAboutAxis(vec3 vec, vec3 axis, float angle){
  // Rotate the vector vec around the given axis by the given angle
   float x = axis.x;
   float y = axis.y;
   float z = axis.z;
  // # S = np.matrix([ [0, -z, y],
  // #                 [z, 0, -x],
  // #                 [-y, x, 0]])
  // # S2 = S @ S
  // #
  // # R = np.identity(3) + np.sin(angle) * S + (1-np.cos(angle)) * S2

  float ca = cos(angle);
  float sa = sin(angle);
  float u = 1 - ca;
  
  float x2  = x*x;
  float y2  = y*y;
  float z2  = z*z;

  // mat3 R = mat3( ca + x2*u,         x*y*u - z*sa,       x*z*u + y*sa
  //                 x*y*u + z*sa,      ca + y2*u,          y*z*u - x*sa,
  //                 x*z*u - y*sa,      y*z*u + x*sa,       ca + z2*u);

  mat3 R = mat3(  ca + x2*u,      x*y*u + z*sa,   x*z*u - y*sa, //1st col
                  x*y*u - z*sa,   ca + y2*u,      y*z*u + x*sa, //2nd col
                  x*z*u + y*sa,   y*z*u - x*sa,   ca + z2*u);   //3rd col

  return R * vec;
}

float GetScatteringAngle(int sca_id){
    return acos(1 - 2 * RandomNumber(vec2(ray_index, sca_id)));
}

float GetPlaneRotationAngle(int sca_id){
  return 2 * PI * RandomNumber(vec2(ray_index, sca_id));
}

vec4 BackwardPropagate(){
    
    bool is_scattered = false;
    vec2 cross_sections;

    int ray_state = VALID;

    int sca_ID = 0;
    float altitude  = history[0].altitude;
    vec3 vec        = history[0].vec;
    float segment_length = history[0].segment_length;
    float sca_angle = history[0].sca_angle;
    float plane_rotation_angle = history[0].plane_rotation_angle;

    vec3 total_vec = vec3(0, 0, 0);

    // while the ray is valid == dont touch the ground or escape in the sky.
    while(ray_state == VALID && sca_ID < scattering_limit){
        // While the ray is not scattered == while the ray goes in a straight line == while we are on the same path segment
        while (ray_state == VALID && !is_scattered){

            altitude += (vec.x * increment_length);  //Increment the altitude
            segment_length += increment_length;      //Increment the segment length

            if (altitude < min_altitude){
                ray_state = TOUCH_GROUND;
                WriteHistory(sca_ID, altitude, vec, -1, -1, -1, segment_length, vec3(0,0,0), ray_state, 0);
            }

            if(altitude > max_altitude){
                ray_state = ESCAPE_IN_SPACE;
                WriteHistory(sca_ID, altitude, vec, -1, -1, -1, segment_length, vec3(0,0,0), ray_state, 0);
            }

            cross_sections = GetAtmosphereCrossSection(altitude) * increment_length;
            float total_cross_section = cross_sections.x + cross_sections.y;

            is_scattered = IsScattered(total_cross_section, sca_ID);
        }

        // Increment the position of the ray with respect to the instrument (origin)
        total_vec += vec * segment_length;

        // Check if the ray is below the max distance limit. If not, return an INVALID ray.
        if (sqrt(total_vec.z*total_vec.z + total_vec.y*total_vec.y) > distance_limit){
          ray_state = INVALID;
        }

        // Compute the scattering angle and change in direction vector of the ray if it has been scattered.
        if (ray_state == VALID){
            sca_angle = GetScatteringAngle(sca_ID);
            plane_rotation_angle = GetPlaneRotationAngle(sca_ID);
            vec3 scattering_plane_normal = GetScatteringPlane(vec, plane_rotation_angle);
            vec3 vec = RotateAboutAxis(vec, scattering_plane_normal, sca_angle);
            //Save the ray history for the given scattering (or end point)
            sca_ID += 1;
            WriteHistory(sca_ID, altitude, vec, sca_angle, cross_sections.x, cross_sections.y, segment_length, scattering_plane_normal, ray_state, plane_rotation_angle);
            segment_length = 0;
        }

    }
    nb_scattering_events = sca_ID;
    return vec4(total_vec, ray_state);

}

mat4 GetScatteringMatrix(float angle){

  float a = 2.118 * 1e-29;  //polarizability in m3. Isotropic.

  float A = a*a; // 5*A and B as defined in de Hulst 1981 p79. isotropic polarizability -> A == B.
  float ca = cos(angle);
  float sa2 = pow(sin(angle), 2);

  float m1 = A - 0.5 * A * sa2;
  float m2 =   - 0.5 * A * sa2;
  float m3 = A * (1 - 0.5 * sa2);
  float m4 = A * ca;

  float Cst = pow(2 * PI * 1e-9 / wavelength, 4);

  mat4 M = mat4(m1, m2, 0,  0,      // 1 col
                m2, m3, 0,  0,      // 2 col
                0,  0,  m4, 0,      // 3 col
                0,  0,  0,  m4);    // 4 col

  return M * Cst;
}

mat4 GetPlaneRotationMatrix(float angle){

  mat4 M = mat4(1, 0,             0,            0,      // 1 col
                0, cos(2*angle), -sin(2*angle), 0,      // 2 col
                0, sin(2*angle), cos(2*angle),  0,      // 3 col
                0, 0,            0,             1);     // 4 col
  return M;
}

float GetPlaneRotationAngle(vec3 P1, vec3 P2){
    float angle = acos(dot(P1, P2));
    return angle;

}

vec4 ForwardPropagate(){
  vec4 stokes_param = vec4(1, 0, 0, 0);

  int ray_state = VALID;
  int sca_id = nb_scattering_events;
  float plane_rotation_angle;
  mat4 L;
  mat4 M;

  while(ray_state == VALID && sca_id >= 0){
    stokes_param.x /= history[sca_id + 1].segment_length * history[sca_id + 1].segment_length;

    plane_rotation_angle = history[sca_id - 1].plane_rotation_angle;
    
    L = GetPlaneRotationMatrix(plane_rotation_angle);
    stokes_param = L * stokes_param;

    M = GetScatteringMatrix(history[sca_id].sca_angle);
    stokes_param = M * stokes_param;

    sca_id -= 1;
    ray_state = history[sca_id].state;

  }

  return stokes_param;
}

////////////////////////////////////////////////////////////////////////////////
// Main function
////////////////////////////////////////////////////////////////////////////////


void main()
{

    uint V_index    = 3 * ray_index;
    uint Vcos_index = V_index + 1;
    uint Vsin_index = V_index + 2;

    uint debug_index = 5 * ray_index;

    int ray_state = INVALID;
    vec4 result;
    vec3 final_position;
    // while (ray_state != TOUCH_GROUND){
        InitializeHistory();
        result = BackwardPropagate();
        ray_state = int(result.w);
        final_position = result.xyz;
        // ray_state = TOUCH_GROUND;
    // }

    //Forward scattering with unity initial intensity condition. Each ray will be scaled later in the python code following the intensity of its starting pos.
    vec4 stokes_param = ForwardPropagate();

    // stokes_param *= wavelength;
    // stokes_param *= distance_limit;


    // Registering output buffers.
    // One mandatory for the stokes paramters of the ray. Similar to the ground and sky scattering shaders. Size 3 * N_rays of floats. Or N_rays of vec4 (with one useless param).
    // If possible, one for the history of each ray to ease the debug, analysis and validation, comparison procedure. Size N_rays * max_sca_events of RayHistoryData if possible.
    observation_data.data[V_index]    = stokes_param.x;
    observation_data.data[Vcos_index] = stokes_param.y;
    observation_data.data[Vsin_index] = stokes_param.z;


    // 5 debug paramters can be saved here for analysis.
    // If you want more debug slots, you must change: 
    //    - "debug_index" at the start of this main function, 
    //    - at the initialization of "debug_data" in simulation.ComputeMultipleScatteringGPU() 
    //    - in the ShaderWrapMS.SaveResults() at the reshape method.
    debug_data.data[debug_index]   = ray_index; 
    debug_data.data[debug_index+1] = ray_state; 
    debug_data.data[debug_index+2] = 0; 
    debug_data.data[debug_index+3] = 0; 
    debug_data.data[debug_index+4] = stokes_param.z; 

}