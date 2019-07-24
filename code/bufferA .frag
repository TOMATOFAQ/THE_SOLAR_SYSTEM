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

const ivec2 txPosL = ivec2(0, 0);
const ivec2 txPosR = ivec2(1, 0);
const ivec2 txPosD = ivec2(2, 0);
const ivec2 txPosU = ivec2(3, 0);
const ivec2 txPospu = ivec2(4, 0);
const ivec2 txPospd = ivec2(5, 0);

// const ivec2 txPosL    = ivec2(4,0);
// const ivec4 txBricks     = ivec4(0,1,13,12);

float LEFT = 0.0;
float RIGHT = 0.0;
float DOWN = 0.0;
float UP = 0.0;
float PUP = 0.0;
float PUD = 0.0;

const float KEY_RIGHT = 39.0;
const float KEY_LEFT = 37.0;
const float KEY_UP = 38.0;
const float KEY_DOWN = 40.0;
const float KPUP = 33.0;
const float KPUD = 34.0;

vec4 loadValue(in ivec2 re) { return texelFetch(iChannel1, re, 0); }

void storeValue(in ivec2 re, in vec4 va, inout vec4 fragColor, in ivec2 p) {
  fragColor = (p == re) ? va : fragColor;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  ivec2 ipx = ivec2(fragCoord - 0.5);

  float PosL = loadValue(txPosL).x;
  float PosR = loadValue(txPosR).x;
  float PosD = loadValue(txPosD).x;
  float PosU = loadValue(txPosU).x;
  float Pospu = loadValue(txPospu).x;
  float Pospd = loadValue(txPospd).x;

  LEFT += texelFetch(iChannel0, ivec2(KEY_LEFT, 0), 0).x + PosL;
  RIGHT += texelFetch(iChannel0, ivec2(KEY_RIGHT, 0), 0).x + PosR;
  DOWN += texelFetch(iChannel0, ivec2(KEY_DOWN, 0), 0).x + PosD;
  UP += texelFetch(iChannel0, ivec2(KEY_UP, 0), 0).x + PosU;
  PUP += texelFetch(iChannel0, ivec2(KPUP, 0), 0).x + Pospu;
  PUD += texelFetch(iChannel0, ivec2(KPUD, 0), 0).x + Pospd;
  fragColor = vec4(0.0);

  storeValue(txPosL, vec4(LEFT, 0.0, 0.0, 0.0), fragColor, ipx);
  storeValue(txPosR, vec4(RIGHT, 0.0, 0.0, 0.0), fragColor, ipx);
  storeValue(txPosD, vec4(DOWN, 0.0, 0.0, 0.0), fragColor, ipx);
  storeValue(txPosU, vec4(UP, 0.0, 0.0, 0.0), fragColor, ipx);
  storeValue(txPospu, vec4(PUP, 0.0, 0.0, 0.0), fragColor, ipx);
  storeValue(txPospd, vec4(PUD, 0.0, 0.0, 0.0), fragColor, ipx);
}
