//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define PROP  (3 * .87f)                      /* \kappa from Regge trajectories */
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
                    __global float4* mom_old,
                    __global float4* pos_new,
                    __global float4* mom_new,
                    __global float4* color,
                    __global float4* cum_force,
                    float dt,
                    uint PARTICLE_NUMBER) 
{
  const unsigned int i = get_global_id(0);

  float4 position = pos_old[i];
  float4 momentum = mom_old[i];
  float mass = momentum.w;
  momentum.w = 0.f;

  
  const float4 c = color[i];


  float4 force = (float4) (0.f, 0.f, 0.f, 0.f);
  for (uint j = 0; j < PARTICLE_NUMBER; ++j) {
    const float4 other_pos = pos_old[j];
    other_pos.w = 0.f;
    const float4 other_col = color[j];
    force += PROP * normalize(position - other_pos) * cornell(c, other_col);
  }

  // force.w should be 0 here - hence momentum.w rests 0
  momentum = momentum + (force * dt * 0.5f);
  
  float4 velocity = momentum / sqrt(dot(momentum, momentum) + mass*mass); 
  velocity.w = 0.f;                     /* this should not be neccessary */

  position += velocity * dt * 0.5f;                          /* half time step to new position */

  /* next half timestep */
  momentum = mom_old[i];
  mass = momentum.w;
  momentum.w = 0.f;
    
  force = (float4) (0.f, 0.f, 0.f, 0.f);
  for (uint j = 0; j < PARTICLE_NUMBER; ++j) {
    const float4 other_pos = pos_old[j];
    other_pos.w = 0.f;
    const float4 other_col = color[j];
    force += PROP * normalize(position - other_pos) * cornell(c, other_col);
  }

  momentum = momentum + (force * dt);   /* in second half step go full dt! */
  velocity = momentum / sqrt(dot(momentum,momentum) + mass*mass); 
  velocity.w = 0.f;


  position = pos_old[i] + velocity * dt;
  
  position.w = 1.f;
  momentum.w = mass;
  pos_new[i] = position;
  mom_new[i] = momentum;
  cum_force[i] = force;
}
