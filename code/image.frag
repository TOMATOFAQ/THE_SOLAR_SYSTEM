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
// CONSTANTS
//============================================================================

//**************************************************
// for track
//**************************************************

#define TrackLength 5
const int NUM_TRIANGLES = TrackLength * 9 * 2;

struct Triangle_t {
  vec3 v0;
  vec3 v1;
  vec3 v2;
  int materialID;
};

Triangle_t Triangle[NUM_TRIANGLES];

//**************************************************
// BACKGROUND
//**************************************************
#define iterations 17
#define formuparam 0.53

#define volsteps 20
#define stepsize 0.1

#define zoom 0.800
#define tile 0.850
#define speed 0.0001

#define brightness 0.0015
#define darkmatter 0.300
#define distfading 0.730
#define saturation 0.850

#define CLIP_PLANE_DIST (1.1 * CAR_TOP_DIST)

//**************************************************
// KEYBOARD CONTROL
//**************************************************

#define PI 3.0

const ivec2 txPosL = ivec2(0, 0);
const ivec2 txPosR = ivec2(1, 0);
const ivec2 txPosD = ivec2(2, 0);
const ivec2 txPosU = ivec2(3, 0);
const ivec2 txPospu = ivec2(4, 0);
const ivec2 txPospd = ivec2(5, 0);

const float EYE_MIN_LATITUDE = -85.0;
const float EYE_MAX_LATITUDE = -85.0;
const float EYE_LATITUDE_INCR = 2.0;
const float EYE_LONGITUDE_INCR = 2.0;
const float EYE_MIN_DIST = 2.0;
const float EYE_INIT_DIST = 20.0;
const float EYE_DIST_INCR = 1.0;

float eyeLatitude = 0.0;
float eyeLongitude = 0.0;
float eyeDistance = EYE_INIT_DIST;

//**************************************************
// planet index
//**************************************************
const int NUM_PLANETS = 10; // 加了多少颗星这里改多少个
const int SUN_INDEX = 0;
const int MERCURY_INDEX = 1;
const int VENUS_INDEX = 2;
const int EARTH_INDEX = 3;
const int MARS_INDEX = 4;
const int JUPITER_INDEX = 5;
const int SATURN_INDEX = 6;
const int URANUS_INDEX = 7;
const int NEPTUNE_INDEX = 8;
const int MOON_INDEX = 9;
const int CORONA_INDEX = 10;
const float SUN_RADIUS = 0.7;
const float CORONA_RADIUS = 0.7;

//**************************************************
// for render planet
//**************************************************
#define saturate(oo) clamp(oo, 0.0, 1.0)

// Quality Settings
int MarchSteps = 16;
// Scene Settings
vec3 ExpPosition = vec3(0.0);
float Radius = 1.0;
vec4 Background = vec4(0.0, 0.0, 0.0, 1.0);
// Noise Settings
const int NoiseSteps = int(4.0);
const float NoiseAmplitude = 0.06;
const float NoiseFrequency = 48.0;
const vec3 Animation = vec3(0.0, -3.0, 0.5);
// Colour Gradient
vec4 BaseColor = vec4(1.0, 1.0, 1.0, 1.0);        // 底色
vec4 InnerFrameColor = vec4(1.0, 0.8, 0.2, 1.0);  // 内焰
vec4 OuterFrameColor = vec4(1.0, 0.03, 0.0, 1.0); // 外焰
vec4 ParticleColor = vec4(0.5, 0.2, 0.2, 1.0);    // 外层颗粒
int PlanetTextureIndex = 0;                       // 选择的材质

//**************************************************
// for ray traycer
//**************************************************
const int NUM_LIGHTS = 1;
const int NUM_MATERIALS = 3;
const int NUM_PLANES = 2;
const int NUM_SPHERES = 0;

const vec3 BACKGROUND_COLOR = vec3(0.0, 0.0, 0.0);
const float FOVY = 50.0 * 3.1415926535 / 180.0;
const float DEFAULT_TMIN = 10.0e-4;
const float DEFAULT_TMAX = 10.0e6;
const int NUM_ITERATIONS = 2;

const float AU = 5.0;
const float V = 1.0;
const float R = 1.0;

//============================================================================
// Define new struct types.
//============================================================================

struct Ray_t {
  vec3 o; // Ray Origin.
  vec3 d; // Ray Direction. A unit vector.
};

struct Sphere_t {
  vec3 center;
  float radius;
  int materialID;
};

struct Light_t {
  vec3 position; // Point light 3D position.
  vec3 I_a;      // For Ambient.
  vec3 I_source; // For Diffuse and Specular.
};

struct Material_t {
  vec3 k_a;  // Ambient coefficient.
  vec3 k_d;  // Diffuse coefficient.
  vec3 k_r;  // Reflected specular coefficient.
  vec3 k_rg; // Global reflection coefficient.
  float n;   // The specular reflection exponent. Ranges from 0.0 to 128.0.
};
struct Planet_t {
  Sphere_t sphere;
  vec4 BaseColor;
  vec4 InnerFrameColor;
  vec4 OuterFrameColor;
  vec4 ParticleColor;
  int PlanetTextureIndex;
};

Planet_t Planet[NUM_PLANETS];
Light_t Light[NUM_LIGHTS];
Material_t Material[NUM_MATERIALS];

//============================================================================
// Functions
//============================================================================

//**************************************************
//行星渲染 START
//**************************************************

void SelectPlanet(int index) {
  ExpPosition = Planet[index].sphere.center;
  Radius = Planet[index].sphere.radius;
  BaseColor = Planet[index].BaseColor;
  InnerFrameColor = Planet[index].InnerFrameColor;
  OuterFrameColor = Planet[index].OuterFrameColor;
  ParticleColor = Planet[index].ParticleColor;
  PlanetTextureIndex = Planet[index].PlanetTextureIndex;
}

float noise(in vec3 x) {
  vec3 p = floor(x);
  vec3 f = fract(x);
  f = f * f * (3.0 - 2.0 * f);
  vec2 uv = (p.xy + vec2(37.0, 17.0) * p.z) + f.xy;
  vec2 rg;
  if (PlanetTextureIndex == 0)
    rg = textureLod(iChannel0, (uv + 0.5) / 256.0, 0.0).yx;
  if (PlanetTextureIndex == 1)
    rg = textureLod(iChannel1, (uv + 0.5) / 256.0, 0.0).yx;

  return -1.0 + 1.7 * mix(rg.x, rg.y, f.z);
}

float Turbulence(vec3 position, float minFreq, float maxFreq, float qWidth) {
  float value = 0.0;
  float cutoff = clamp(0.5 / qWidth, 0.0, maxFreq);
  float fade;
  float fOut = minFreq;
  for (int i = NoiseSteps; i >= 0; i--) {
    if (fOut >= 0.5 * cutoff)
      break;
    fOut *= 2.0;
    value += abs(noise(position * fOut)) / fOut;
  }
  fade = clamp(2.0 * (cutoff - fOut) / cutoff, 0.0, 1.0);
  value += fade * abs(noise(position * fOut)) / fOut;
  return 1.0 - value;
}

float SphereDist(vec3 position) {
  return length(position - ExpPosition) - Radius;
}

vec4 Shade(float distance) {
  float c1 = saturate(distance * 5.0 + 0.5);
  float c2 = saturate(distance * 5.0);
  float c3 = saturate(distance * 3.4 - 0.5);
  vec4 a = mix(BaseColor, InnerFrameColor, c1);
  vec4 b = mix(a, OuterFrameColor, c2);
  return mix(b, ParticleColor, c3);
}

// Draws the scene
float RenderScene(vec3 position, out float distance) {
  float noise = Turbulence(position * NoiseFrequency + Animation * iTime * 0.24,
                           0.1, 1.5, 0.03) *
                NoiseAmplitude;
  noise = saturate(abs(noise));
  distance = SphereDist(position) - noise;
  return noise;
}

// Basic ray marching method.
vec4 March(vec3 rayOrigin, vec3 rayStep,
           out vec3 position) { // 根据射线原点和射线方向返回颜色
  position = rayOrigin;
  float distance;
  float displacement;
  for (int step = MarchSteps; step >= 0; --step) {
    displacement = RenderScene(
        position,
        distance); // 根据当前位置进行渲染，获得一个 float 类型的噪音值
    if (distance <
        0.05) // 如果在规定步长内，离噪声的距离已经小于0.05，那么break
      break;
    position +=
        rayStep *
        distance; // ❓，方向乘距离前进是什么操作，这种前进方式决定的是什么？
  }
  return mix(
      Shade(displacement), Background,
      float(distance >= 0.5)); // 一个开关，根据距离噪声的距离返回是背景还是颜色
}

bool IntersectSphere(vec3 ro, vec3 rd, vec3 pos, float radius, out float t,
                     out vec3 intersectPoint, out vec3 hitNormal,
                     out vec3 color) // 行星渲染中检测是否撞到粒子
{
  vec3 relDistance = (ro - pos);
  float b = dot(relDistance, rd);
  float c = dot(relDistance, relDistance) - radius * radius;
  float d = b * b - c;
  if (d <= 0.0)
    return false;
  else {
    color = vec3(March(ro, rd, intersectPoint)); // 最后居然是这里错了？？
    t = length(intersectPoint - ro);
    hitNormal = normalize(intersectPoint - pos);
    intersectPoint = ro + rd * (-b - sqrt(d));
    // intersectPoint返回的是所打中粒子最终的点

    return true;
  }
}

//**************************************************
//轨迹渲染 START
//**************************************************

vec3 getPlanetLocation(int planet_index, float time) {
  if (planet_index == SUN_INDEX)
    return vec3(0.0, 0.0, 0.0);
  if (planet_index == MERCURY_INDEX)
    return vec3(0.55 * AU * cos(time + iTime * 1.8 * V),
                0.55 * AU * sin(time + iTime * 1.8 * V), 0.0);
  if (planet_index == VENUS_INDEX)
    return vec3(0.73 * AU * cos(time + iTime * 1.45 * V),
                0.73 * AU * sin(time + iTime * 1.45 * V), 0.0);
  if (planet_index == EARTH_INDEX)
    return vec3(AU * cos(time + iTime * V), AU * sin(time + iTime * V), 0.0);
  if (planet_index == MARS_INDEX)
    return vec3(1.34 * AU * cos(time + iTime * 0.87 * V),
                1.34 * AU * sin(time + iTime * 0.87 * V), 0.0);
  if (planet_index == JUPITER_INDEX)
    return vec3(1.78 * AU * cos(time + iTime * 0.66 * V),
                1.78 * AU * sin(time + iTime * 0.66 * V), 0.0);
  if (planet_index == SATURN_INDEX)
    return vec3(2.34 * AU * cos(time + iTime * 0.45 * V),
                2.34 * AU * sin(time + iTime * 0.45 * V), 0.0);
  if (planet_index == URANUS_INDEX)
    return vec3(2.89 * AU * cos(time + iTime * 0.3 * V),
                2.89 * AU * sin(time + iTime * 0.3 * V), 0.0);
  if (planet_index == NEPTUNE_INDEX)
    return vec3(3.24 * AU * cos(time + iTime * 0.2 * V),
                3.24 * AU * sin(time + iTime * 0.2 * V), 0.0);
  if (planet_index == MOON_INDEX)
    return vec3(
        (0.50 * R * cos(time + iTime * 12.0 * V) + AU * cos(time + iTime * V)),
        (0.50 * R * sin(time + iTime * 12.0 * V) + AU * sin(time + iTime * V)),
        0.0);
}

void generatePlanetTrack(int planet_index) {
  // the first triangle of this track was position in Triangle[planet_index *
  // TrackLength * 2]
  int loc = int(planet_index * TrackLength * 2); // 节省性能
  float radius = Planet[planet_index].sphere.radius * 1.5;
  // vec3 currentPosition = getPlanetLocation(planet_index, 0.0);
  vec3 currentPosition = getPlanetLocation(planet_index, 0.0);
  vec3 nextPosition;
  for (int i = 0; i <= TrackLength; i++) {

    // Don't ask me why, it just work
    float m = float(i);
    float n = m + 1.0;
    nextPosition = getPlanetLocation(planet_index, -float(n) * 0.5); // i 秒前
    float o = n + 1.0;
    if ((i % 2) == 1) {
      Triangle[loc].v1 = currentPosition + vec3(0, radius / n, 4.0 * m);
      Triangle[loc].v2 = currentPosition + vec3(0, -radius / n, 4.0 * m);
      Triangle[loc].v0 = nextPosition + vec3(radius / o, 0, 4.0 * n);
      Triangle[loc].materialID = 0;
      loc++;
      Triangle[loc].v1 = nextPosition + vec3(radius / o, 0, 4.0 * n);
      Triangle[loc].v2 = nextPosition + vec3(-radius / o, 0, 4.0 * n);
      Triangle[loc].v0 = currentPosition + vec3(0, radius / n, 4.0 * m);
      Triangle[loc].materialID = 0;
      loc++;
    } else {
      Triangle[loc].v1 = currentPosition + vec3(radius / n, 0, 4.0 * m);
      Triangle[loc].v2 = currentPosition + vec3(-radius / n, 0, 4.0 * m);
      Triangle[loc].v0 = nextPosition + vec3(0, radius / o, 4.0 * n);
      Triangle[loc].materialID = 0;
      loc++;
      Triangle[loc].v1 = nextPosition + vec3(0, radius / o, 4.0 * n);
      Triangle[loc].v2 = nextPosition + vec3(0, -radius / o, 4.0 * n);
      Triangle[loc].v0 = currentPosition + vec3(radius / n, 0, 4.0 * m);
      Triangle[loc].materialID = 0;
      loc++;
    }
    currentPosition = nextPosition;
  }
}

//**************************************************
// InitScene START
//**************************************************

void InitScene() {
  // 不反射的 material.
  Material[0].k_d = vec3(0.0);
  Material[0].k_a = vec3(0.0);
  Material[0].k_r = vec3(0.0);
  Material[0].k_rg = vec3(0.0);
  Material[0].n = 0.0;

  // 略微反射 material.
  Material[1].k_d = vec3(0.5, 0.5, 0.5);
  Material[1].k_a = 0.2 * Material[0].k_d;
  Material[1].k_r = 2.0 * Material[0].k_d;
  Material[1].k_rg = 0.5 * Material[0].k_r;
  Material[1].n = 64.0;

  // Light 0.
  Light[0].position = vec3(0.0, 0.0, 0.0);
  Light[0].I_a = vec3(1.0, 1.0, 1.0);
  Light[0].I_source = vec3(1.0, 1.0, 1.0);

  Planet[SUN_INDEX].sphere.center = vec3(0.0, 0.0, 0.0); // 以此控制运动
  Planet[SUN_INDEX].sphere.radius = 2.0 * R;
  Planet[SUN_INDEX].sphere.materialID = 0;
  Planet[SUN_INDEX].BaseColor = vec4(1.0, 1.0, 1.0, 1.0);
  Planet[SUN_INDEX].InnerFrameColor = vec4(1.0, 0.8, 0.2, 1.0);
  Planet[SUN_INDEX].OuterFrameColor = vec4(1.0, 0.03, 0.0, 1.0);
  Planet[SUN_INDEX].ParticleColor == vec4(0.5, 0.2, 0.2, 1.0);
  Planet[SUN_INDEX].PlanetTextureIndex = 0; // RGBA noise Medium

  Planet[MERCURY_INDEX].sphere.center = vec3(
      0.55 * AU * cos(iTime * 1.8 * V), 0.55 * AU * sin(iTime * 1.8 * V), 0.0);
  Planet[MERCURY_INDEX].sphere.radius = 0.08 * R;
  Planet[MERCURY_INDEX].sphere.materialID = 1;
  Planet[MERCURY_INDEX].BaseColor = vec4(0.95686, 0.82353, 0.38039, 1.0);
  Planet[MERCURY_INDEX].InnerFrameColor = vec4(1.0, 0.8, 0.2, 1.0);
  Planet[MERCURY_INDEX].OuterFrameColor = vec4(0.62, 0.58, 0.50, 1.0);
  Planet[MERCURY_INDEX].ParticleColor == vec4(0.62, 0.58, 0.50, 1.0);
  Planet[MERCURY_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[VENUS_INDEX].sphere.center =
      vec3(0.73 * AU * cos(iTime * 1.45 * V), 0.73 * AU * sin(iTime * 1.45 * V),
           0.0);
  Planet[VENUS_INDEX].sphere.radius = 0.24 * R;
  Planet[VENUS_INDEX].sphere.materialID = 1;
  Planet[VENUS_INDEX].BaseColor = vec4(0.61961, 0.79216, 0.37255, 1.0);
  Planet[VENUS_INDEX].InnerFrameColor = vec4(1.0, 0.8, 0.2, 1.0);
  Planet[VENUS_INDEX].OuterFrameColor = vec4(1.0, 0.89, 0.51, 1.0);
  Planet[VENUS_INDEX].ParticleColor == vec4(0.92, 0.52, 0.12, 1.0);
  Planet[VENUS_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[EARTH_INDEX].sphere.center =
      vec3(AU * cos(iTime * V), AU * sin(iTime * V), 0.0);
  Planet[EARTH_INDEX].sphere.radius = 0.26 * R;
  Planet[EARTH_INDEX].sphere.materialID = 1;
  Planet[EARTH_INDEX].BaseColor = vec4(0.88235, 0.81569, 0.61961, 1.0);
  Planet[EARTH_INDEX].InnerFrameColor = vec4(0.03, 0.18, 0.32, 1.0);
  Planet[EARTH_INDEX].OuterFrameColor = vec4(0.03, 0.18, 0.54, 1.0);
  Planet[EARTH_INDEX].ParticleColor == vec4(1.0, 1.0, 1.0, 1.0);
  Planet[EARTH_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[MARS_INDEX].sphere.center =
      vec3(1.34 * AU * cos(iTime * 0.87 * V), 1.34 * AU * sin(iTime * 0.87 * V),
           0.0);
  Planet[MARS_INDEX].sphere.radius = 0.15 * R;
  Planet[MARS_INDEX].sphere.materialID = 1;
  Planet[MARS_INDEX].BaseColor = vec4(0.88235, 0.81569, 0.61961, 1.0);
  Planet[MARS_INDEX].InnerFrameColor = vec4(1.0, 0.38, 0.0, 1.0);
  Planet[MARS_INDEX].OuterFrameColor = vec4(1.0, 0.38, 0.0, 1.0);
  Planet[MARS_INDEX].ParticleColor == vec4(0.36, 0.14, 0.70, 1.0);
  Planet[MARS_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[JUPITER_INDEX].sphere.center =
      vec3(1.78 * AU * cos(iTime * 0.66 * V), 1.78 * AU * sin(iTime * 0.66 * V),
           0.0);
  Planet[JUPITER_INDEX].sphere.radius = 0.5 * R;
  Planet[JUPITER_INDEX].sphere.materialID = 1;
  Planet[JUPITER_INDEX].BaseColor = vec4(0.83922, 0.67059, 0.53333, 1.0);
  Planet[JUPITER_INDEX].InnerFrameColor = vec4(0.96, 0.87, 0.70, 1.0);
  Planet[JUPITER_INDEX].OuterFrameColor = vec4(0.96, 0.87, 0.70, 1.0);
  Planet[JUPITER_INDEX].ParticleColor == vec4(0.98, 0.90, 0.78, 1.0);
  Planet[JUPITER_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[SATURN_INDEX].sphere.center =
      vec3(2.34 * AU * cos(iTime * 0.45 * V), 2.34 * AU * sin(iTime * 0.45 * V),
           0.0);
  Planet[SATURN_INDEX].sphere.radius = 0.46 * R;
  Planet[SATURN_INDEX].sphere.materialID = 1;
  Planet[SATURN_INDEX].BaseColor = vec4(0.89412, 0.87059, 0.63922, 1.0);
  Planet[SATURN_INDEX].InnerFrameColor = vec4(0.49, 1.0, 0.0, 1.0);
  Planet[SATURN_INDEX].OuterFrameColor = vec4(0.49, 1.0, 0.0, 1.0);
  Planet[SATURN_INDEX].ParticleColor == vec4(0.23, 0.56, 0.25, 1.0);
  Planet[SATURN_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[URANUS_INDEX].sphere.center = vec3(
      2.89 * AU * cos(iTime * 0.3 * V), 2.89 * AU * sin(iTime * 0.3 * V), 0.0);
  Planet[URANUS_INDEX].sphere.radius = 0.22 * R;
  Planet[URANUS_INDEX].sphere.materialID = 1;
  Planet[URANUS_INDEX].BaseColor = vec4(0.00784, 0.87059, 0.95686, 1.0);
  Planet[URANUS_INDEX].InnerFrameColor = vec4(0.25, 0.95, 0.81, 1.0);
  Planet[URANUS_INDEX].OuterFrameColor = vec4(0.25, 0.95, 0.81, 1.0);
  Planet[URANUS_INDEX].ParticleColor == vec4(0.25, 0.95, 0.81, 1.0);
  Planet[URANUS_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[NEPTUNE_INDEX].sphere.center = vec3(
      3.24 * AU * cos(iTime * 0.2 * V), 3.24 * AU * sin(iTime * 0.2 * V), 0.0);
  Planet[NEPTUNE_INDEX].sphere.radius = 0.27 * R;
  Planet[NEPTUNE_INDEX].sphere.materialID = 1;
  Planet[NEPTUNE_INDEX].BaseColor = vec4(0.00000, 0.51765, 0.89020, 1.0);
  Planet[NEPTUNE_INDEX].InnerFrameColor = vec4(0.25, 0.41, 0.88, 1.0);
  Planet[NEPTUNE_INDEX].OuterFrameColor = vec4(0.25, 0.41, 0.88, 1.0);
  Planet[NEPTUNE_INDEX].ParticleColor == vec4(0.25, 0.41, 0.88, 1.0);
  Planet[NEPTUNE_INDEX].PlanetTextureIndex = 1; // Rusty Metal

  Planet[MOON_INDEX].sphere.center =
      vec3((0.50 * R * cos(iTime * 12.0 * V) + AU * cos(iTime * V)),
           (0.50 * R * sin(iTime * 12.0 * V) + AU * sin(iTime * V)), 0.0);
  Planet[MOON_INDEX].sphere.radius = 0.008 * R;
  Planet[MOON_INDEX].sphere.materialID = 1;
  Planet[MOON_INDEX].BaseColor = vec4(0.00000, 0.51765, 0.89020, 1.0);
  Planet[MOON_INDEX].InnerFrameColor = vec4(0.96, 0.87, 0.70, 1.0);
  Planet[MOON_INDEX].OuterFrameColor = vec4(0.96, 0.87, 0.70, 1.0);
  Planet[MOON_INDEX].ParticleColor == vec4(0.98, 0.90, 0.78, 1.0);
  Planet[MOON_INDEX].PlanetTextureIndex = 1; // Rusty Metal

 // generatePlanetTrack(SUN_INDEX);
  generatePlanetTrack(MERCURY_INDEX);
  generatePlanetTrack(VENUS_INDEX);
  generatePlanetTrack(EARTH_INDEX);
  generatePlanetTrack(MARS_INDEX);
  generatePlanetTrack(JUPITER_INDEX);
  generatePlanetTrack(SATURN_INDEX);
  generatePlanetTrack(URANUS_INDEX);
  generatePlanetTrack(NEPTUNE_INDEX);
  // generatePlanetTrack();
}

bool IntersectTriangle(in Triangle_t triangle, in Ray_t ray, in float tmin,
                       in float tmax, out float t, out vec3 hitPos,
                       out vec3 hitNormal) {

  vec3 A = triangle.v0;
  vec3 B = triangle.v1;
  vec3 C = triangle.v2;

  mat3 mat_of_beta =
      mat3(A.x - ray.o.x, A.y - ray.o.y, A.z - ray.o.z, A.x - C.x, A.y - C.y,
           A.z - C.z, ray.d.x, ray.d.y, ray.d.z);

  mat3 mat_of_gama =
      mat3(A.x - B.x, A.y - B.y, A.z - B.z, A.x - ray.o.x, A.y - ray.o.y,
           A.z - ray.o.z, ray.d.x, ray.d.y, ray.d.z);

  mat3 solve_det = mat3(A.x - B.x, A.y - B.y, A.z - B.z, A.x - C.x, A.y - C.y,
                        A.z - C.z, ray.d.x, ray.d.y, ray.d.z);

  float beta = determinant(mat_of_beta) / determinant(solve_det);
  float gama = determinant(mat_of_gama) / determinant(solve_det);

  if (beta > 0.0 && gama > 0.0 && beta + gama < 1.0) {
    float temp_t;
    mat3 mat_of_t =
        mat3(A.x - B.x, A.y - B.y, A.z - B.z, A.x - C.x, A.y - C.y, A.z - C.z,
             A.x - ray.o.x, A.y - ray.o.y, A.z - ray.o.z);
    temp_t = determinant(mat_of_t) / determinant(solve_det);
    if (temp_t > tmin && temp_t < tmax) {
      t = temp_t;
      hitPos = ray.o + t * ray.d;
      hitNormal = normalize(cross(B - A, C - A));
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

bool IntersectSphere(in Sphere_t sph, in Ray_t ray, in float tmin,
                     in float tmax, out float t, out vec3 hitPos,
                     out vec3 hitNormal) {
  vec3 rayLocalPosition = ray.o - sph.center;
  float a = dot(ray.d, ray.d);
  float b = 2.0 * dot(ray.d, rayLocalPosition);
  float c = dot(rayLocalPosition, rayLocalPosition) - pow(sph.radius, 2.0);
  float d = pow(b, 2.0) - 4.0 * a * c;
  float t1 = (-b + sqrt(d)) / (2.0 * a);
  float t2 = (-b - sqrt(d)) / (2.0 * a);
  if (d < 0.0)
    return false;

  t = tmax;
  if (t1 > 0.0 && t1 < t)
    t = t1;
  if (t2 > 0.0 && t2 < t)
    t = t2;
  if (t <= tmin || t >= tmax)
    return false;

  hitPos = ray.o + t * ray.d;
  hitNormal = normalize(hitPos - sph.center);

  return true;
}

//**************************************************
// Ray tracing START
//**************************************************

vec3 PhongLighting(in vec3 L, in vec3 N, in vec3 V, in bool inShadow,
                   in Material_t mat, in Light_t light) {

  if (inShadow) {
    return light.I_a * mat.k_a;
  } else {
    vec3 R = reflect(-L, N);
    float N_dot_L = max(0.0, dot(N, L));
    float R_dot_V = max(0.0, dot(R, V));
    float R_dot_V_pow_n = (R_dot_V == 0.0) ? 0.0 : pow(R_dot_V, mat.n);

    return light.I_a * mat.k_a +
           light.I_source * (mat.k_d * N_dot_L + mat.k_r * R_dot_V_pow_n);
  }
}

vec3 CastRay(in Ray_t ray, out bool hasHit, out vec3 hitPos, out vec3 hitNormal,
             out vec3 k_rg) {
  // Find whether and where the ray hits some object.
  // Take the nearest hit point.
  vec3 I_local = vec3(0.0); // Result color will be accumulated here.
  vec3 I_temp;

  bool hasHitSomething = false;
  float nearest_t =
      DEFAULT_TMAX;       // The ray parameter t at the nearest hit point.
  vec3 nearest_hitPos;    // 3D position of the nearest hit point.
  vec3 nearest_hitNormal; // Normal vector at the nearest hit point.
  int nearest_hitMatID;   // MaterialID of the object at the nearest hit point.

  float temp_t;
  vec3 temp_hitPos;
  vec3 temp_hitNormal;
  bool temp_hasHit;

  // 检测是否射到行星，若有顺便求其颜色
  for (int i = 0; i < NUM_PLANETS; i++) {
    SelectPlanet(i);
    if (IntersectSphere(ray.o, ray.d, Planet[i].sphere.center,
                        Planet[i].sphere.radius + NoiseAmplitude * 12.0, temp_t,
                        temp_hitPos, temp_hitNormal, I_temp)) {
      if (I_temp.r > 0.0 || I_temp.g > 0.0) {
        // 判定击中的是否太阳的实体。可能确实是太阳的球但是是为了占位而设置成的空
        if (!hasHitSomething || temp_t < nearest_t) {

          nearest_t = temp_t;
          nearest_hitPos = temp_hitPos;
          nearest_hitNormal = temp_hitNormal;
          nearest_hitMatID = Planet[i].sphere.materialID;
          hasHitSomething = true;

          I_local = I_temp;
        }
      }
    }
  }

  // 单纯在上面加一层纱
  for (int i = 0; i < NUM_TRIANGLES; i++) {
    if (IntersectTriangle(Triangle[i], ray, DEFAULT_TMIN, DEFAULT_TMAX, temp_t,
                          temp_hitPos, temp_hitNormal)) {
      if (!hasHitSomething || temp_t < nearest_t) {
        hasHitSomething = true;
        I_local += vec3(0.1);
        break;
      }
    }
  }

  // One of the output results.
  hasHit = hasHitSomething;
  if (!hasHitSomething)
    return BACKGROUND_COLOR;

  for (int i = 0; i < NUM_LIGHTS; i++) {
    Ray_t shadowRay;
    shadowRay.o = Light[i].position;
    shadowRay.d = normalize(nearest_hitPos - Light[i].position);

    bool inShadow = false;

    for (int j = 0; j < NUM_PLANETS; j++) { // 阴影检测
      if (j == SUN_INDEX)
        continue;
      if (IntersectSphere(Planet[j].sphere, shadowRay, DEFAULT_TMIN,
                          DEFAULT_TMAX, temp_t, temp_hitPos, temp_hitNormal))
        if (length(nearest_hitPos - Light[j].position) - DEFAULT_TMIN >
            length(temp_hitPos - Light[j].position)) {
          inShadow = true;
          break;
        }
    }

    I_local += PhongLighting(-shadowRay.d, nearest_hitNormal,
                             normalize(ray.o - nearest_hitPos), inShadow,
                             Material[nearest_hitMatID], Light[i]);
  }

  // Populate output results.
  hitPos = nearest_hitPos;
  hitNormal = nearest_hitNormal;
  k_rg = Material[nearest_hitMatID].k_rg;

  return I_local;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  InitScene();

  float LEFT = texelFetch(iChannel3, txPosL, 0).x;
  float RIGHT = texelFetch(iChannel3, txPosR, 0).x;
  float DOWN = texelFetch(iChannel3, txPosD, 0).x;
  float UP = texelFetch(iChannel3, txPosU, 0).x;
  float PUP = texelFetch(iChannel3, txPospu, 0).x;
  float PUD = texelFetch(iChannel3, txPospd, 0).x;

  // BACKGROUND
  vec2 uv = fragCoord.xy / iResolution.xy - .5;
  uv.y *= iResolution.y / iResolution.x;
  vec3 dir = vec3(uv * zoom, 1.);
  float time = iTime * speed + .25;

  // mouse rotation
  float a1 = .5 + iMouse.x / iResolution.x * 2.;
  float a2 = .8 + iMouse.y / iResolution.y * 2.;
  mat2 rot1 = mat2(cos(a1), sin(a1), -sin(a1), cos(a1));
  mat2 rot2 = mat2(cos(a2), sin(a2), -sin(a2), cos(a2));
  dir.xz *= rot1;
  dir.xy *= rot2;
  vec3 from = vec3(1.0, 0.5, 0.5);
  from += vec3(time * 2., time, -2.);
  from.xz *= rot1;
  from.xy *= rot2;

  // volumetric rendering
  float s = 0.1, fade = 1.;
  vec3 v = vec3(0.);
  for (int r = 0; r < volsteps; r++) {
    vec3 p = from + s * dir * .5;
    //       vec3 p=s*dir*.5;
    p = abs(vec3(tile) - mod(p, vec3(tile * 2.))); // tiling fold
    float pa, a = pa = 0.;
    for (int i = 0; i < iterations; i++) {
      p = abs(p) / dot(p, p) - formuparam; // the magic formula
      a += abs(length(p) - pa);            // absolute sum of average change
      pa = length(p);
    }
    float dm = max(0., darkmatter - a * a * .001); // dark matter
    a *= a * a;                                    // add contrast
    if (r > 6)
      fade *= 1. - dm; // dark matter, don't render near
    // v+=vec3(dm,dm*.5,0.);
    v += fade;
    v += vec3(s, s * s, s * s * s * s) * a * brightness *
         fade;          // coloring based on distance
    fade *= distfading; // distance fading
    s += stepsize;
  }
  v = mix(vec3(length(v)), v, saturation); // color adjust
  fragColor = vec4(v * 0.01, 1.);

  // Scale pixel 2D position such that its y coordinate is in [-1.0, 1.0].
  vec2 pixel_pos = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.y; //

  // Position the camera.
  //   vec3 cam_pos = vec3(0.0,0.0,20.0);
  vec3 cam_pos = vec3(
      (eyeDistance + (PUP - PUD) * EYE_DIST_INCR) *
          sin((eyeLongitude + (RIGHT - LEFT) * EYE_LONGITUDE_INCR) * PI /
              180.0) *
          cos((eyeLatitude + (DOWN - UP) * EYE_LATITUDE_INCR) * PI / 180.0),
      (eyeDistance + (PUP - PUD) * EYE_DIST_INCR) *
          sin((eyeLatitude + (DOWN - UP) * EYE_LATITUDE_INCR) * PI / 180.0),
      (eyeDistance + (PUP - PUD) * EYE_DIST_INCR) *
          cos((eyeLongitude + (RIGHT - LEFT) * EYE_LONGITUDE_INCR) * PI /
              180.0) *
          cos((eyeLatitude + (DOWN - UP) * EYE_LATITUDE_INCR) * PI / 180.0));

	float temp_distance = 1.0;
    if(iTime >= 20.0 && iTime <=25.0)
    {
    	temp_distance =(1.0-(iTime - 20.0)/5.0);
		cam_pos = cam_pos * temp_distance;
    }
  vec3 cam_lookat = vec3(0.0, 0.0, 0.0);
  vec3 cam_up_vec = vec3(0.0, 1.0, 0.0);

  // Set up camera coordinate frame in world space.
  // 视角坐标系的基在世界坐标系下的表达
  vec3 cam_z_axis = normalize(cam_pos - cam_lookat);
  vec3 cam_x_axis = normalize(cross(cam_up_vec, cam_z_axis));
  vec3 cam_y_axis = normalize(cross(cam_z_axis, cam_x_axis));

  // Create primary ray.
  float pixel_pos_z = -1.0 / tan(FOVY / 2.0);
  Ray_t pRay;
  pRay.o = cam_pos;
  pRay.d = normalize(pixel_pos.x * cam_x_axis + pixel_pos.y * cam_y_axis +
                     pixel_pos_z * cam_z_axis);

  vec3 I_result = vec3(0.0);
  vec3 compounded_k_rg = vec3(1.0); // 混合的全局反射系数
  Ray_t nextRay = pRay;
  Ray_t ray = pRay;

  for (int level = 0; level <= NUM_ITERATIONS; level++) {
    bool hasHit;
    vec3 hitPos, hitNormal, k_rg;
    vec3 I_local = CastRay(nextRay, hasHit, hitPos, hitNormal, k_rg);
    I_result += compounded_k_rg * I_local;
    if (!hasHit)
      break;
    compounded_k_rg *= k_rg; // 复合的光照变量
    nextRay = Ray_t(hitPos, normalize(reflect(nextRay.d, hitNormal)));
  }
  fragColor += vec4(I_result, 1.0);
    if(iTime >25.0)
  {
  	vec2 uv = fragCoord.xy / iResolution.xy;
	fragColor = texture(iChannel2, uv);
  }
}
