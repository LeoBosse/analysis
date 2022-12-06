#version 430

// Set up our compute groups
layout(local_size_x=COMPUTE_SIZE_X, local_size_y=COMPUTE_SIZE_Y) in;

// Input uniforms go here if you need them.
// Some examples:
uniform float instrument_azimut;
uniform float instrument_elevation;
uniform vec3  instrument_los;
uniform float instrument_area;
uniform float instrument_fov;
uniform float wavelength;
uniform float los_nb_points;
uniform float los_altitudes;
uniform float atm_altitudes;
uniform float atm_abs;
uniform float map_azimuts;
uniform float map_distances;


// Structure of the ball data
// struct Ball
// {
//     vec4 pos;
//     vec4 vel;
//     vec4 color;
// };

// Input buffer
layout(std430, binding=0) buffer emissions
{
    vec3 emissions[];
} Emissions;

// Output buffer
layout(std430, binding=1) buffer measures
{
    vec3 measures[];
} Measures;


// Compute the index of the emission point in the azimut list, elevation list and scattering point list.
vec3 GetAzElLOSFromIndex(int shader_index){
    int nb_az = map_azimuts.length();
    int nb_dist = map_distances.length();

    int az_index = int(shader_index / (nb_dist * los_nb_points));

    int dist_index = int((shader_index - az_index * nb_dist * los_nb_points) / los_nb_points);

    int los_index = shader_index % los_nb_points;

    return vec3(az_index, dist_index, los_index);
}




void main()
{
    int shader_index = int(gl_GlobalInvocationID);
    vec3 indexes = GetAzElLOSFromIndex(shader_index);
    int iaz = indexes.x;
    int idist = indexes.y;
    int ilos = indexes.z;
    float azimut = map_azimuts[iaz];
    float distance = map_distances[idist];
    float sca_altitude = los_altitudes[ilos];


    vec3 emission = Emissions.emissions[shader_index];

    Measures.measures[shader_index] = emission;
}
