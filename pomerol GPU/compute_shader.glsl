#version 440



#define X %%local_sizeX%%
#define Y %%local_sizeY%%
#define Z %%local_sizeZ%%


// Set up our compute groups
layout(local_size_x = X, local_size_y = Y, local_size_z = Z) in;

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform float instrument_altitude;
uniform float instrument_lon;
uniform float instrument_lat;
uniform float instrument_area;
uniform float instrument_fov;
uniform float wavelength;
uniform int   los_nb_points;
uniform float map_delta_az;
uniform bool  is_point_source;
uniform int   atm_nb_altitudes;
uniform int   atm_nb_angles;
uniform bool  use_aerosol;

uniform float test[100];


// Structure of the ball data
// struct Ball
// {
//     vec4 pos;
//     vec4 vel;
//     vec4 color;
// };


// Stokes paramter buffers
layout(std430, binding=0) buffer in_0{
    float in_V[];};
layout(std430, binding=1) buffer in_1{
    float in_Vcos[];};
layout(std430, binding=2) buffer in_2{
    float in_Vsin[];};

layout(std430, binding=3) buffer out_3{
    float out_V[];};
layout(std430, binding=4) buffer out_4{
    float out_Vcos[];};
layout(std430, binding=5) buffer out_5{
    float out_Vsin[];};



// Compute the index of the emission point in the azimut list, elevation list and scattering point list.
/* vec3 GetAzElLOSFromIndex(int shader_index){
    int nb_az = map_azimuts.length();
    int nb_dist = map_distances.length();

    int az_index = int(shader_index / (nb_dist * los_nb_points));

    int dist_index = int((shader_index - az_index * nb_dist * los_nb_points) / los_nb_points);

    int los_index = shader_index % los_nb_points;

    return vec3(az_index, dist_index, los_index);
} */




void main()
{
    const int buffer_index = int(gl_GlobalInvocationID.x * gl_NumWorkGroups.y + gl_GlobalInvocationID.y);
    /*
    vec3 indexes = GetAzElLOSFromIndex(shader_index);
    int iaz = indexes.x;
    int idist = indexes.y;
    int ilos = indexes.z;
    float azimut = map_azimuts[iaz];
    float distance = map_distances[idist];
    float sca_altitude = los_altitudes[ilos];


    vec3 emission = Emissions.emissions[shader_index];
*/
    out_V[buffer_index] = test[buffer_index];
}
