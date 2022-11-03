#ifndef LTC_INPUT_INCLUDED
#define LTC_INPUT_INCLUDED

TEXTURE2D(_BaseMap);            SAMPLER(sampler_BaseMap);
TEXTURE2D(_BumpMap);     
TEXTURE2D(_RoughnessMap);
TEXTURE2D(_LTC_MatrixTexture);
TEXTURE2D(_LTC_MagnitueTexture);
TEXTURE2D_ARRAY(_FilteredLightTexture);
StructuredBuffer<float3> _PolygonalLightVertexPos;
CBUFFER_START(UnityPerMaterial)
float4 _BaseMap_ST;
half4 _BaseColor;
half4 _SpecularColor;
float _Roughness;
float _Intensity;
float _SpecularIntensity;
float  _IsTwoSided = false;
int _Layer = 0;
const float LUT_SIZE = 64.0;
const float LUT_SCALE = (64.0 - 1.0) / 64.0;
const float LUT_BIAS = 0.5 / 64.0;
const float _PI = 3.14159265;
CBUFFER_END


#endif