//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define PROP (3*0.87f)                      /* \kappa from Regge trajectories */
#define EPS 1e-5f

/*----------------------------------------------------------------------
  color is (r,g,b,alpha)
  anticolor is bit-flipped, i.e. r = (1,0,0,1), -r = (0,1,1,1)
  i.e. transparency remains unchanged
  ---------------------------------------------------------------------*/

constant float4 unit = (float4)(1.f, 1.f, 1.f, 0.f); //ignore alpha channel 

inline float cornell(float4 my_color, float4 its_color)
{
  const float4 combination = its_color + my_color;

  if ( fast_length(my_color - its_color) <= EPS )
    return +1.f/3.f;                    /* q + q */
  if ( dot(unit, combination) == 2.f )
    return -1.f/6.f;                    /* r+g; r+b; g+b */
  if ( dot(unit, combination) == 4.f )
    return -1.f/6.f;                    /* -r-g; -r-b; -g-b */
  if ( (combination.x == 1.f) && (combination.y == 1.f) && (combination.z == 1.f) )
    return -1.f/3.f;                    /* q + -q */
  return 1.f/6.f;                       /* -rg; -rb; b-g; b-r; g-r; g-b */

}  

__kernel void nbody(__global float4* pos_old, 
                    __global float4* vel_old,
                    __global float4* pos_new,
                    __global float4* vel_new,
                    __global float4* color, 
                    float dt,
                    uint PARTICLE_NUMBER) 
{
  const unsigned int i = get_global_id(0);

  float4 p = pos_old[i];
  float mean = 0.f;

  for (uint j = 0; j < PARTICLE_NUMBER; ++j)
    mean += fast_distance(p, pos_old[j]);

  mean /= PARTICLE_NUMBER;


  float4 v = vel_old[i];
  float mass = v.w;
  v.w = 0.f;
  p.w = 0.f;
  
  const float4 c = color[i];


  float4 force = (float4) (0.f, 0.f, 0.f, 0.f);
  float4 p_neu = (float4) (0.f, 0.f, 0.f, 0.f);

  for (uint j = 0; j < PARTICLE_NUMBER; ++j) {
    const float4 other_pos = pos_old[j];
    other_pos.w = 0.f;
    const float4 other_col = color[j];
    force += PROP * normalize(p - other_pos) * ( cornell(c, other_col) );//+ log(length(p-other_pos)+EPS) );
  }

  p_neu = (mass * v * 1.f/sqrt(EPS+1.f - length(v)*length(v)))+(force * dt * 0.5f);

  v = p_neu / sqrt(length(p_neu)*length(p_neu) + mass*mass); 
  v.w = 0.f;

  p += v * dt * 0.5f;                          /* half time step to new position */

  /* next half timestep */
  v = vel_old[i];
  mass = v.w;
  v.w = 0.f;
    
  force = (float4) (0.f, 0.f, 0.f, 0.f);
  for (uint j = 0; j < PARTICLE_NUMBER; ++j) {
    const float4 other_pos = pos_old[j];
    other_pos.w = 0.f;
    const float4 other_col = color[j];
    force += PROP * normalize(p - other_pos) * ( cornell(c, other_col));// + log(length(p-other_pos)+EPS) );
  }
  p_neu = (mass * v * 1.f/sqrt(EPS+1.f - length(v)*length(v)))+(force * dt);
  
  v = p_neu / sqrt(length(p_neu)*length(p_neu) + mass*mass); 
  v.w = 0.f;


  p = pos_old[i] + v * dt;
  
  p.w = 1.f;
  v.w = mass;
  pos_new[i] = p;
  vel_new[i] = v;
}
