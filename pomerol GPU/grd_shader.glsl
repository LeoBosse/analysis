#version 440

// Input uniforms go here if you need them.
// Some examples:

uniform float instrument_azimut;
uniform float instrument_elevation;

uniform sampler2DRect instrument_los;
uniform float instrument_area;
uniform float instrument_fov;
uniform float wavelength;
uniform int   los_nb_points;
uniform sampler2DRect los_altitudes;
uniform sampler2DRect atm_altitudes;
uniform sampler2DRect atm_abs;

uniform sampler2DRect map_azimuts;
uniform sampler2DRect map_distances;



#define X %%local_sizeX%%
#define Y %%local_sizeY%%
#define Z %%local_sizeZ%%

// Set up our compute groups
layout(local_size_x = X, local_size_y = Y, local_size_z = Z) in;


// Structure of the ball data
// struct Ball
// {
//     vec4 pos;
//     vec4 vel;
//     vec4 color;
// };

// Input buffer
layout(std430, binding=0) buffer in_0{
    float in_V[];};
layout(std430, binding=1) buffer in_1{
    float in_Vcos[];};
layout(std430, binding=2) buffer in_2{
    float in_Vsin[];};

// Output buffer
layout(std430, binding=3) buffer out_3{
    float out_V[];};
layout(std430, binding=4) buffer out_4{
    float out_Vcos[];};
layout(std430, binding=5) buffer out_5{
    float out_Vsin[];};


float ReadTexture(sampler2DRect text, int index){
  return texture(text, vec2(index, 0.)).x;
}

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
    int buffer_index = int(gl_GlobalInvocationID.x * gl_NumWorkGroups.y + gl_GlobalInvocationID.y);


    /*vec3 indexes = GetAzElLOSFromIndex(shader_index);

    const int iaz = int(indexes.x);
    const int idist = int(indexes.y);
    const int ilos = int(indexes.z); */

    /* float azimut = map_azimuts[iaz];
    float distance = map_distances[idist];
    float sca_altitude = los_altitudes[ilos]; */

    out_V[buffer_index] = ReadTexture(los_altitudes, buffer_index);
    /* out_Vcos[buffer_index] = in_Vcos[buffer_index];
    out_Vsin[buffer_index] = in_Vsin[buffer_index]; */
}
