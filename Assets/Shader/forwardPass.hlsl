#ifndef FORWARD_PASS_INCLUDED
#define FORWARD_PASS_INCLUDED
#include "ltc.hlsl"
struct Attributes
{
    float4 positionOS   : POSITION;
    float3 normalOS     : NORMAL;
    float4 tangentOS    : TANGENT;
    float2 texcoord     : TEXCOORD0;
};

struct Varyings
{
    float4 positionCS               : SV_POSITION;
    float2 uv : TEXCOORD0;
    float3 positionWS               : TEXCOORD1;
    float4 normalWS                 : TEXCOORD2;    // xyz: normal, w: viewDir.x
    float4 tangentWS                : TEXCOORD3;    // xyz: tangent, w: viewDir.y
    float4 bitangentWS              : TEXCOORD4;
    half4 fogFactorAndVertexLight   : TEXCOORD5;
    float4 shadowCoord              : TEXCOORD6;
    float3 viewDirWS              : TEXCOORD7;
};

Varyings ForwardPassVertex(Attributes input)
{
    Varyings output = (Varyings)0;

    VertexPositionInputs vertexInput = GetVertexPositionInputs(input.positionOS.xyz);
    VertexNormalInputs normalInput = GetVertexNormalInputs(input.normalOS, input.tangentOS);
    float3 viewDirWS = GetCameraPositionWS() - vertexInput.positionWS;
    output.normalWS = float4(normalInput.normalWS, viewDirWS.x);
    output.tangentWS = float4(normalInput.tangentWS, viewDirWS.y);
    output.bitangentWS = float4(normalInput.bitangentWS, viewDirWS.z);
    //output.normalWS.xyz =  TransformObjectToWorldNormal(input.normalOS.xyz,true);
    //output.tangentWS.xyz = TransformObjectToWorldDir(input.tangentOS);
    //output.bitangentWS.xyz = cross(output.normalWS.xyz,output.tangentWS.xyz) * input.tangentOS.w * unity_WorldTransformParams.w;
    
    float3 vertexLight = VertexLighting(vertexInput.positionWS, normalInput.normalWS);
    half fogFactor = ComputeFogFactor(vertexInput.positionCS.z);
    output.uv = TRANSFORM_TEX(input.texcoord, _BaseMap);
    output.viewDirWS = viewDirWS;

    OUTPUT_LIGHTMAP_UV(input.lightmapUV, unity_LightmapST, output.lightmapUV);
    output.fogFactorAndVertexLight = half4(fogFactor, vertexLight);

    output.positionWS = vertexInput.positionWS;

    output.shadowCoord = GetShadowCoord(vertexInput);
    output.positionCS = vertexInput.positionCS;

    return output;
}

half4 ForwardPassFragment(Varyings input) : SV_Target
{

    half4 baseColor = SAMPLE_TEXTURE2D(_BaseMap, sampler_BaseMap,input.uv);
    float3 n = UnpackNormal(SAMPLE_TEXTURE2D(_BumpMap, sampler_BaseMap,input.uv));
    n = TransformTangentToWorld(n,real3x3(input.tangentWS.xyz, input.bitangentWS.xyz, input.normalWS.xyz));
    //float3 n = normalize(input.normalWS.xyz);	//其他几何体的话应该由其法线乘以模型矩阵来算
    //float3 n = float3(0.0,1.0,0.0);
	float3 v = normalize(GetCameraPositionWS() - input.positionWS);
	float r = SAMPLE_TEXTURE2D(_RoughnessMap, sampler_BaseMap,input.uv).x;
	float3 posWorld = input.positionWS;

	float2 UV = LTC_Coords(dot(n, v),max( _Roughness + r,0.0));
	
	//test
	UV = float2(UV.x ,1-UV.y);
    
	float4 LTCMatrixComponents = SAMPLE_TEXTURE2D(_LTC_MatrixTexture, sampler_BaseMap,UV);
	float3x3 LTCMatrix = float3x3
	(
		float3(1, 0, LTCMatrixComponents.y),
		float3(0, LTCMatrixComponents.z, 0),
		float3(LTCMatrixComponents.w, 0, LTCMatrixComponents.x)
	);
	LTCMatrix = transpose(LTCMatrix);
	
    float3x3 DiffuseLTCMatrix = float3x3
	(
		float3(1.0, 0.0, 0.0),
		float3(0.0, 1.0, 0.0),
		float3(0.0, 0.0, 1.0)
	);
	float3 Diffuse = integrateLTC(n, v, posWorld, DiffuseLTCMatrix, _PolygonalLightVertexPos, _IsTwoSided);
	float3 Specular = integrateLTC(n, v, posWorld, LTCMatrix, _PolygonalLightVertexPos, _IsTwoSided);
	float w = SAMPLE_TEXTURE2D(_LTC_MagnitueTexture,sampler_BaseMap, UV).w;
	Specular *= float3(w,w,w);
	float2 Schlick = SAMPLE_TEXTURE2D(_LTC_MagnitueTexture,sampler_BaseMap, UV).xy;
	Specular *= _SpecularColor.xyz * Schlick.x + (float3(1.0,1.0,1.0) - _SpecularColor.xyz) * Schlick.y;

	float3 ResultColor = _Intensity * (baseColor.xyz * Diffuse*_BaseColor.xyz + Specular*_SpecularIntensity);
	//float3 ResultColor = u_Intensity * Diffuse * u_DiffuseColor;
	ResultColor /=  PI;
//	ResultColor /= 2.0 * PI;
	half4 color = float4(ResultColor, 1.0);
	//half4 color = float4(baseColor.xyz * Diffuse, 1.0);
    return color;
    //return half4(_PolygonalLightVertexPos[3].xyz,1);
    //return half4(baseColor.xyz * SAMPLE_TEXTURE2D_ARRAY(_FilteredLightTexture ,sampler_BaseMap, input.uv ,0.0).rgb ,1.0);
    //return half4(baseColor.xyz * SAMPLE_TEXTURE2D(_LTC_MatrixTexture ,sampler_BaseMap, input.uv).rgb ,1.0);
}
#endif