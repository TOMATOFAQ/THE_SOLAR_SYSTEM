//============================================================================

// GROUP NUMBER:13

//

// STUDENT NAME:WU PENGYU

// NUS User ID.: t0918566

//

// STUDENT NAME: WANG SIYUAN

// NUS User ID.: t0918594

//

// STUDENT NAME: FAN QIANYI	

// NUS User ID.: t0918683
//============================================================================
// Constants.
//============================================================================
const int NUM_LIGHTS = 3;
const int NUM_MATERIALS = 7;
const int NUM_PLANES = 4;
const int NUM_SPHERES = 2;
const int NUM_TRIANGLE = 8;

const vec3 BACKGROUND_COLOR = vec3( 0.1, 0.2, 0.6 );

// Vertical field-of-view angle of camera. In radians.
const float FOVY = 50.0 * 3.1415926535 / 180.0; 

// Use this for avoiding the "epsilon problem" or the shadow acne problem. // �Ա���������
const float DEFAULT_TMIN = 10.0e-4;

// Use this for tmax for non-shadow ray intersection test.
const float DEFAULT_TMAX = 10.0e6;

// Equivalent to number of recursion levels (0 means ray-casting only). 
// We are using iterations to replace recursions.
// ����׷�ٵ�������
const int NUM_ITERATIONS = 3;


//============================================================================
// Define new struct types.
//============================================================================
struct Ray_t {
    vec3 o;  // Ray Origin.
    vec3 d;  // Ray Direction. A unit vector.
};

struct Plane_t {
    // The plane equation is Ax + By + Cz + D = 0.
    float A, B, C, D;
    int materialID;
};

struct Sphere_t {
    vec3 center;
    float radius;
    int materialID;
};

struct Triangle_t {
    vec3 v0;
	vec3 v1;
	vec3 v2;
    int materialID;
};

struct Light_t {
    vec3 position;  // Point light 3D position.
    vec3 I_a;       // For Ambient.
    vec3 I_source;  // For Diffuse and Specular.
};

struct Material_t {
    vec3 k_a;   // Ambient coefficient.
    vec3 k_d;   // Diffuse coefficient.
    vec3 k_r;   // Reflected specular coefficient.
    vec3 k_rg;  // Global reflection coefficient. // ȫ�ֹ��գ���ɶ�õģ����ڵ�������⣿ I = I_local  +  k_rg * I_reflected
    float n;    // The specular reflection exponent. Ranges from 0.0 to 128.0. 
};


//============================================================================
// Global scene data.
//============================================================================
Plane_t Plane[NUM_PLANES];
Sphere_t Sphere[NUM_SPHERES];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];
Triangle_t Triangle[NUM_TRIANGLE];

/////////////////////////////////////////////////////////////////////////////
// Initializes the scene.
/////////////////////////////////////////////////////////////////////////////
void InitScene()
{
    // Horizontal plane.
    Plane[0].A = 0.0;
    Plane[0].B = 1.0;
    Plane[0].C = 0.0;
    Plane[0].D = 0.0;
    Plane[0].materialID = 5;

    // Vertical plane.
    Plane[1].A = 0.0;
    Plane[1].B = 0.0;
    Plane[1].C = 1.0;
    Plane[1].D = 3.5;
    Plane[1].materialID = 4;
    
    //left floor
    Plane[2].A = 2.8;
    Plane[2].B = 0.0;
    Plane[2].C = 1.0;
    Plane[2].D = 12.0;
    Plane[2].materialID = 4;
    
	//up floor
    Plane[3].A = 0.0;
    Plane[3].B = 1.0;
    Plane[3].C = 0.0;
    Plane[3].D = -5.0;
    Plane[3].materialID = 6;    
    
    // Center bouncing sphere.
    Sphere[0].center = vec3( 3.0 , abs(sin(0.3 * iTime)) + 0.7+1.0, 4.0 ); // �Դ˿����˶� abs(sin(2.0 * iTime)) + 0.7+1.0
    Sphere[0].radius = 1.0;
    Sphere[0].materialID = 1;

    // Circling sphere.
    Sphere[1].center = vec3( 3.0 + 2.5  * cos(iTime), 1.7+ abs(sin(0.3 * iTime)) , 4.0 + 2.5 * sin(iTime) );
    Sphere[1].radius = 0.5;
    Sphere[1].materialID = 2;

	//Triangle 
	Triangle[0].v0 = vec3(4.0, 0.0, 4.0);
	Triangle[0].v1 = vec3(3.0, 1.0, 4.0);
	Triangle[0].v2 = vec3(3.0, 0.0, 5.0);
	Triangle[0].materialID = 1;
    
    Triangle[1].v0 = vec3(4.0, 0.0, 4.0);
	Triangle[1].v1 = vec3(3.0, 1.0, 4.0);
	Triangle[1].v2 = vec3(3.0, 0.0, 5.0);
	Triangle[1].materialID = 1;

    Triangle[2].v0 = vec3(3.0, 0.0, 3.0);
	Triangle[2].v1 = vec3(3.0, 1.0, 4.0);
	Triangle[2].v2 = vec3(4.0, 0.0, 4.0);
	Triangle[2].materialID = 1;

    Triangle[3].v0 = vec3(2.0, 0.0, 4.0);
	Triangle[3].v1 = vec3(3.0, 1.0, 4.0);
	Triangle[3].v2 = vec3(3.0, 0.0, 3.0);
	Triangle[3].materialID = 1;
    
    Triangle[4].v0 = vec3(2.0, 5.0, 4.0);
	Triangle[4].v1 = vec3(3.0 +cos(2.0*iTime) , 6.0, 3.0 + sin(2.0 * iTime));
	Triangle[4].v2 = vec3(4.0, 4.0, 4.0);
	Triangle[4].materialID = 5;
    
    Triangle[5].v0 = vec3(2.0, 3.0, 4.0);
	Triangle[5].v1 = vec3(1.5 + 2.0*cos(iTime), 5.0, 1.0+sin(iTime));
	Triangle[5].v2 = vec3(3.0, 1.0, 4.0);
	Triangle[5].materialID = 5;
    
    Triangle[6].v0 = vec3(5.5, 4.5+ sin(iTime), 7.0+ cos(iTime));
	Triangle[6].v1 = vec3(5.0, 6.0, 7.0);
	Triangle[6].v2 = vec3(6.0 + sin(iTime), 4.0, 6.0 +cos(iTime));
	Triangle[6].materialID = 6;
    
    Triangle[7].v0 = vec3(7.0 + 2.0*sin(iTime), 3.0+ cos(iTime), 3.0);
	Triangle[7].v1 = vec3(8.0, 3.5+ cos(iTime), 3.0 + sin(iTime));
	Triangle[7].v2 = vec3(9.0+ cos(iTime), 2.5 +cos(iTime), 2.0);
	Triangle[7].materialID = 5;

    // Silver material.
    Material[0].k_d = vec3( 0.5, 0.5, 0.5 );
    Material[0].k_a = 0.2 * Material[0].k_d;
    Material[0].k_r = 2.0 * Material[0].k_d;
    Material[0].k_rg = 0.5 * Material[0].k_r;
    Material[0].n = 64.0;

    // Gold material.
    Material[1].k_d = vec3( 0.8, 0.7, 0.3 );
    Material[1].k_a = 0.2 * Material[1].k_d;
    Material[1].k_r = 2.0 * Material[1].k_d;
    Material[1].k_rg = 0.5 * Material[1].k_r;
    Material[1].n = 64.0;

    // Green plastic material.
    Material[2].k_d = vec3( 0.0, 0.8, 0.0 );
    Material[2].k_a = 0.2 * Material[2].k_d;
    Material[2].k_r = vec3( 1.0, 1.0, 1.0 );
    Material[2].k_rg = 0.5 * Material[2].k_r;
    Material[2].n = 128.0;
    
    // mirror material: blue.
    Material[3].k_d = vec3( 0.3, 0.4, 0.4 );
    Material[3].k_a = 1.5 * Material[3].k_d;
    Material[3].k_r = 3.0 * Material[3].k_d;
    Material[3].k_rg = 0.5 * Material[3].k_r;
    Material[3].n = 64.0;
    
	// mirror material:
    Material[4].k_d = vec3( 0.3, 0.3, 0.2 );
    Material[4].k_a = 2. * Material[3].k_d;
    Material[4].k_r = 3.0 * Material[3].k_d;
    Material[4].k_rg = 0.5 * Material[3].k_r;
    Material[4].n = 64.0;
    
    // mirror white material.
    Material[5].k_d = vec3( 1.0, 1.0, 1.0 );
    Material[5].k_a = 2. * Material[3].k_d;
    Material[5].k_r = 3.0 * Material[3].k_d;
    Material[5].k_rg = 0.5 * Material[3].k_r;
    Material[5].n = 64.0;    
    
    // mirror silver material.
    Material[6].k_d = vec3( 0.9, 0.9, 0.9 );
    Material[6].k_a = 2. * Material[3].k_d;
    Material[6].k_r = 3.0 * Material[3].k_d;
    Material[6].k_rg = 0.5 * Material[3].k_r;
    Material[6].n = 64.0;
    


    // Light 0.
    Light[0].position = vec3( -4.0, 3.0, 4.0 );
    Light[0].I_a = vec3( 0.1, 0.1, 0.07 );
    Light[0].I_source = vec3( 0.4, 0.4, 0.4 );

    // Light 1.
    Light[1].position = vec3( 4.5, 2.0, 0.0 );
    Light[1].I_a = vec3( 0.1, 0.1, 0.07 );
    Light[1].I_source = vec3( 0.4, 0.4, 0.4 );
    
    // Light 2.
    Light[2].position = vec3( 10, 2, 8.0 );
    Light[2].I_a = vec3( 0.05, 0.05, 0.08 );
    Light[2].I_source = vec3( 0.25, 0.25, 0.25 );
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is such an intersection, outputs the value of t, the position
// of the intersection (hitPos) and the normal vector at the intersection 
// (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray, in float tmin, in float tmax,
                     out float t, out vec3 hitPos, out vec3 hitNormal ) 
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;

    // We have a hit -- output results.
    t = t0;
    hitPos = ray.o + t0 * ray.d;
    hitNormal = normalize( N );
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a plane and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectPlane( in Plane_t pln, in Ray_t ray, in float tmin, in float tmax )
{
    vec3 N = vec3( pln.A, pln.B, pln.C );
    float NRd = dot( N, ray.d );
    float NRo = dot( N, ray.o );
    float t0 = (-pln.D - NRo) / NRd;
    if ( t0 < tmin || t0 > tmax ) return false;
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
// If there is one or two such intersections, outputs the value of the 
// smaller t, the position of the intersection (hitPos) and the normal 
// vector at the intersection (hitNormal).
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray, in float tmin, in float tmax,
                      out float t, out vec3 hitPos, out vec3 hitNormal ) 
{
    /////////////////////////////////
    // TASK: WRITE YOUR CODE HERE. //
    /////////////////////////////////
    vec3 rayLocalPosition = ray.o-sph.center;
    float a = dot(ray.d,ray.d);
    float b = 2.0*dot(ray.d,rayLocalPosition);
    float c = dot(rayLocalPosition,rayLocalPosition) - pow(sph.radius,2.0);
    float d = pow(b,2.0)-4.0*a*c;
    float t1 = (-b+sqrt(d))/(2.0*a);
    float t2 = (-b-sqrt(d))/(2.0*a);
    if ( d < 0.0 ) return false;
    
    t = tmax;
    if(t1>0.0&&t1<t)
        t = t1;
    if(t2>0.0&&t2<t)
        t = t2;
    if(t<=tmin||t>=tmax)
        return false;
        
    hitPos = ray.o + t*ray.d;
    hitNormal = normalize(hitPos-sph.center);
    
    return true;
}



/////////////////////////////////////////////////////////////////////////////
// Computes intersection between a sphere and a ray.
// Returns true if there is an intersection where the ray parameter t is
// between tmin and tmax, otherwise returns false.
/////////////////////////////////////////////////////////////////////////////
bool IntersectSphere( in Sphere_t sph, in Ray_t ray, in float tmin, in float tmax )
{
    /////////////////////////////////
    // TASK: WRITE YOUR CODE HERE. //
    /////////////////////////////////
    vec3 rayLocalPosition = ray.o-sph.center;
    float a = dot(ray.d,ray.d);
    float b = 2.0*dot(ray.d,rayLocalPosition);
    float c = dot(rayLocalPosition,rayLocalPosition) - pow(sph.radius,2.0);
    float d = pow(b,2.0)-4.0*a*c;
    float t1 = (-b+sqrt(d))/(2.0*a);
    float t2 = (-b-sqrt(d))/(2.0*a);
    if( d < 0.0) return false;
    return true;
}


////////////////////////////////////////////////////////////////////////////
//compute ray tracing in Triangle
////////////////////////////////////////////////////////////////////////////
bool IntersectTriangle( in Triangle_t tri,in Ray_t ray, in float tmin, in float tmax,
						out float t, out vec3 hitPos, out vec3 hitNormal) 
{ 
	float beta, gamma, temp_t;

	float Ax = (tri.v0).x;
	float Ay = (tri.v0).y;
	float Az = (tri.v0).z;

	float Bx = (tri.v1).x;
	float By = (tri.v1).y;
	float Bz = (tri.v1).z;

	float Cx = (tri.v2).x;
	float Cy = (tri.v2).y;
	float Cz = (tri.v2).z;

	float Rox = (ray.o).x;
	float Roy = (ray.o).y;
	float Roz = (ray.o).z;

	float Rdx = (ray.d).x;
	float Rdy = (ray.d).y;
	float Rdz = (ray.d).z;

    mat3 solve_beta = mat3(
	Ax-Rox, Ax-Cx, Rdx,
	Ay-Roy, Ay-Cy, Rdy,
        Az-Roz, Az-Cz, Rdz);

	mat3 solve_gamma = mat3(
	Ax-Bx, Ax-Rox, Rdx,
	Ay-By, Ay-Roy, Rdy,
	Az-Bz, Az-Roz, Rdz);

	mat3 solve_t = mat3(
	Ax-Bx, Ax-Cx, Ax-Rox,
	Ay-By, Ay-Cy, Ay-Roy,
	Az-Bz, Az-Cz, Az-Roz);

	mat3 solve_det = mat3(
	Ax-Bx, Ay-By, Az-Bz,
	Ax-Cx, Ay-Cy, Az-Cz,
	Rdx, Rdy, Rdz);

	beta = determinant(solve_beta)/determinant(solve_det);
	gamma = determinant(solve_gamma)/determinant(solve_det);
	temp_t = determinant(solve_t)/determinant(solve_det);

	if(beta>0.0 && gamma>0.0 && beta+gamma <1.0)
	{
		if(temp_t > tmin && temp_t < tmax)
		{
			hitNormal = normalize(cross(tri.v1-tri.v0, tri.v2-tri.v0));
			hitPos = ray.o + temp_t * ray.d;
			t = temp_t;
			return true;
		}
		else 
		{
			return false;
		}
	}
	else
		return false;
} 

bool IntersectTriangle( in Triangle_t tri, in Ray_t ray, in float tmin, in float tmax) 
{ 
	float beta, gamma, temp_t;

	float Ax = (tri.v0).x;
	float Ay = (tri.v0).y;
	float Az = (tri.v0).z;

	float Bx = (tri.v1).x;
	float By = (tri.v1).y;
	float Bz = (tri.v1).z;

	float Cx = (tri.v2).x;
	float Cy = (tri.v2).y;
	float Cz = (tri.v2).z;

	float Rox = (ray.o).x;
	float Roy = (ray.o).y;
	float Roz = (ray.o).z;

	float Rdx = (ray.d).x;
	float Rdy = (ray.d).y;
	float Rdz = (ray.d).z;

	mat3 solve_beta = mat3(
	Ax-Rox, Ax-Cx, Rdx,
	Ay-Roy, Ay-Cy, Rdy,
	Az-Roz, Az-Cz, Rdz);

	mat3 solve_gamma = mat3(
	Ax-Bx, Ax-Rox, Rdx,
	Ay-By, Ay-Roy, Rdy,
	Az-Bz, Az-Roz, Rdz);

	mat3 solve_t = mat3(
	Ax-Bx, Ax-Cx, Ax-Rox,
	Ay-By, Ay-Cy, Ay-Roy,
	Az-Bz, Az-Cz, Az-Roz);

	mat3 solve_det = mat3(
	Ax-Bx, Ay-By, Az-Bz,
	Ax-Cx, Ay-Cy, Az-Cz,
	Rdx, Rdy, Rdz);

	beta = determinant(solve_beta)/determinant(solve_det);
	gamma = determinant(solve_gamma)/determinant(solve_det);
	temp_t = determinant(solve_t)/determinant(solve_det);

	if(beta>0.0 && gamma>0.0 && beta+gamma <1.0)
	{
		if(temp_t > tmin && temp_t < tmax)
		{
			return true;
		}
		else 
		{
			return false;
		}
	}
	else
		return false;
} 

/////////////////////////////////////////////////////////////////////////////
// Computes (I_a * k_a) + k_shadow * I_source * [ k_d * (N.L) + k_r * (R.V)^n ].
// Input vectors L, N and V are pointing AWAY from surface point.
// Assume all vectors L, N and V are unit vectors.
/////////////////////////////////////////////////////////////////////////////
vec3 PhongLighting( in vec3 L, in vec3 N, in vec3 V, in bool inShadow, 
                    in Material_t mat, in Light_t light )
{
    if ( inShadow ) {
        return light.I_a * mat.k_a;
    }
    else {
        vec3 R = reflect( -L, N );
        float N_dot_L = max( 0.0, dot( N, L ) );
        float R_dot_V = max( 0.0, dot( R, V ) );
        float R_dot_V_pow_n = ( R_dot_V == 0.0 )? 0.0 : pow( R_dot_V, mat.n );

        return light.I_a * mat.k_a + 
               light.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
    }
}


/////////////////////////////////////////////////////////////////////////////
// Casts a ray into the scene and returns color computed at the nearest
// intersection point. The color is the sum of light from all light sources,
// each computed using Phong Lighting Model, with consideration of
// whether the interesection point is being shadowed from the light.
// If there is no interesection, returns the background color, and outputs
// hasHit as false.
// If there is intersection, returns the computed color, and outputs
// hasHit as true, the 3D position of the intersection (hitPos), the
// normal vector at the intersection (hitNormal), and the k_rg value
// of the material of the intersected object.
/////////////////////////////////////////////////////////////////////////////

vec3 CastRay( in Ray_t ray, 
              out bool hasHit, out vec3 hitPos, out vec3 hitNormal, out vec3 k_rg ) 
{
    // Find whether and where the ray hits some object. 
    // Take the nearest hit point.

    bool hasHitSomething = false;
    float nearest_t = DEFAULT_TMAX;   // The ray parameter t at the nearest hit point.
    vec3 nearest_hitPos;              // 3D position of the nearest hit point.
    vec3 nearest_hitNormal;           // Normal vector at the nearest hit point.
    int nearest_hitMatID;             // MaterialID of the object at the nearest hit point.

    float temp_t;
    vec3 temp_hitPos;
    vec3 temp_hitNormal;
    bool temp_hasHit;

    /////////////////////////////////////////////////////////////////////////////
    // TASK:
    // * Try interesecting input ray with all the planes and spheres,
    //   and record the front-most (nearest) interesection.
    // * If there is interesection, need to record hasHitSomething,
    //   nearest_t, nearest_hitPos, nearest_hitNormal, nearest_hitMatID.
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////
    // TASK: WRITE YOUR CODE HERE. //
    /////////////////////////////////

    for(int i =0;i<NUM_PLANES;i++){
        if(IntersectPlane( Plane[i],ray,  DEFAULT_TMIN, DEFAULT_TMAX,temp_t, temp_hitPos,  temp_hitNormal)){
            if(!hasHitSomething || temp_t < nearest_t){
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Plane[i].materialID;
                hasHitSomething = true;
            }
        }
    }

    for(int i=0;i<NUM_SPHERES;i++){
        if(IntersectSphere( Sphere[i],ray,  DEFAULT_TMIN, DEFAULT_TMAX,temp_t, temp_hitPos,  temp_hitNormal)){
            if(!hasHitSomething || temp_t < nearest_t){
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Sphere[i].materialID;
                hasHitSomething = true;
            }
        }
    }

	for(int i=0;i<NUM_TRIANGLE;i++){
        if(IntersectTriangle( Triangle[i], ray, DEFAULT_TMIN, DEFAULT_TMAX, temp_t,temp_hitPos,temp_hitNormal)){
            if(!hasHitSomething || temp_t < nearest_t){
                nearest_t = temp_t;
                nearest_hitPos = temp_hitPos;
                nearest_hitNormal = temp_hitNormal;
                nearest_hitMatID = Triangle[i].materialID;
                hasHitSomething = true;
            }
        }
    }

    // One of the output results.
    hasHit = hasHitSomething;
    if ( !hasHitSomething ) return BACKGROUND_COLOR;

    vec3 I_local = vec3( 0.0 );  // Result color will be accumulated here.

    /////////////////////////////////////////////////////////////////////////////
    // TASK:
    // * Accumulate lighting from each light source on the nearest hit point. 
    //   They are all accumulated into I_local.
    // * For each light source, make a shadow ray, and check if the shadow ray
    //   intersects any of the objects (the planes and spheres) between the 
    //   nearest hit point and the light source.
    // * Then, call PhongLighting() to compute lighting for this light source.
    /////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////
    // TASK: WRITE YOUR CODE HERE. //
    /////////////////////////////////

    for(int i=0;i<NUM_LIGHTS;i++){
        Ray_t shadowRay;
        shadowRay.o = Light[i].position;
        shadowRay.d = normalize(nearest_hitPos - Light[i].position);

        bool inShadow = false;

        for(int j=0;j<NUM_PLANES;j++){
            if(IntersectPlane(Plane[j],shadowRay,DEFAULT_TMIN,DEFAULT_TMAX,temp_t,temp_hitPos,temp_hitNormal))
                if(length(nearest_hitPos-Light[i].position)-DEFAULT_TMIN>length(temp_hitPos-Light[i].position)){
                    inShadow = true;
                    break;
                }
        }
        for(int j=0;j<NUM_SPHERES;j++){
            if(IntersectSphere(Sphere[j],shadowRay,DEFAULT_TMIN,DEFAULT_TMAX,temp_t,temp_hitPos,temp_hitNormal))
                if(length(nearest_hitPos-Light[i].position)-DEFAULT_TMIN>length(temp_hitPos-Light[i].position)){
                    inShadow = true;
                    break;
                }
        }
		
		for(int j=0;j<NUM_TRIANGLE;j++)
        {
            if(IntersectTriangle(Triangle[j], shadowRay, DEFAULT_TMIN, DEFAULT_TMAX))
                if(length(nearest_hitPos-Light[i].position)-DEFAULT_TMIN>length(temp_hitPos-Light[i].position)){
                    inShadow = true;
                    break;
                }
        }
        
        I_local += PhongLighting( -shadowRay.d, nearest_hitNormal,normalize(ray.o - nearest_hitPos), inShadow, 
                    Material[nearest_hitMatID] , Light[i] );
    }

    // Populate output results.
    hitPos = nearest_hitPos;
    hitNormal = nearest_hitNormal;
    k_rg = Material[nearest_hitMatID].k_rg;

    return I_local;
}



/////////////////////////////////////////////////////////////////////////////
// Execution of fragment shader starts here.
// 1. Initializes the scene.
// 2. Compute a primary ray for the current pixel (fragment).
// 3. Trace ray into the scene with NUM_ITERATIONS recursion levels.
/////////////////////////////////////////////////////////////////////////////
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    InitScene();
    
    // Scale pixel 2D position such that its y coordinate is in [-1.0, 1.0].
    vec2 pixel_pos = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.y; 

    // Position the camera. // ���ﻹ����������ϵ
    vec3 cam_pos = vec3( 10.0, 3.0, 10.0 ); 
    vec3 cam_lookat = vec3( 0.0, 2.0, 0.0 );
    vec3 cam_up_vec = vec3( 0.0, 1.0, 0.0 );

    // Set up camera coordinate frame in world space. �ӽ�����ϵ�Ļ�����������ϵ�µı��
    vec3 cam_z_axis = normalize( cam_pos - cam_lookat );
    vec3 cam_x_axis = normalize( cross(cam_up_vec, cam_z_axis));
    vec3 cam_y_axis = normalize( cross(cam_z_axis, cam_x_axis));  
    
    // Create primary ray.
    float pixel_pos_z = -1.0 / tan(FOVY / 2.0); 

    Ray_t pRay;
    pRay.o = cam_pos;
    pRay.d = normalize( pixel_pos.x * cam_x_axis  +  pixel_pos.y * cam_y_axis  +  pixel_pos_z * cam_z_axis );
	
    // Start Ray Tracing.
    // Use iterations to emulate the recursion.

    vec3 I_result = vec3( 0.0 );
    vec3 compounded_k_rg = vec3( 1.0 ); // ��ϵ�ȫ�ַ���ϵ��
    Ray_t nextRay = pRay;

    for ( int level = 0; level <= NUM_ITERATIONS; level++ ) 
    {
        bool hasHit;
        vec3 hitPos, hitNormal, k_rg;

        vec3 I_local = CastRay( nextRay, hasHit, hitPos, hitNormal, k_rg );

        I_result += compounded_k_rg * I_local;

        if ( !hasHit ) break;

        compounded_k_rg *= k_rg; // ���ϵĹ��ձ���

        nextRay = Ray_t( hitPos, normalize( reflect(nextRay.d, hitNormal) ) );
    }

    fragColor = vec4( I_result, 1.0 );
    if(iTime <= 10.0 )
        fragColor = vec4(0.0, 0.0, 0.0, 1.0);
    if(iTime >25.0 && iTime <30.0)
    {
    	fragColor =fragColor + vec4( 1.0-(iTime-25.0)/5.0, 1.0-(iTime-25.0)/5.0, 1.0-(iTime-25.0)/5.0, 0.0 );
    }
}
