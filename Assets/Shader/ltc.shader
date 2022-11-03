Shader "Unlit/ltc"
{
    Properties
    {
        _BaseMap("Texture", 2D) = "white" {}
        _BaseColor("基础色", Color) = (1,1,1,1)
        _SpecularColor("高光色", Color) = (1,1,1,1)
        _BumpMap("法线图", 2D) = "white" {}
        _RoughnessMap("粗糙度图", 2D) = "white" {}
        _Roughness("粗糙度偏移值", Range(-1.0,1.0)) = 0
        _Intensity("强度", Range(0.0,1.0)) = 1
        _SpecularIntensity("高光强度",Range(0.0,1.0)) = 1
        [Toggle]_IsTwoSided("双面", Float) = 0
        _Layer("层数", Int) = 1
        _LTC_MatrixTexture("_LTC_MatrixTexture", 2D) = "black" {}
        _LTC_MagnitueTexture("_LTC_MagnitueTexture", 2D) = "black" {}
        _FilteredLightTexture("法线图", 2DArray) = "white" {}
    }
    SubShader
    {
        Tags
        {
            "RenderPipeline"="UniversalPipeline"
            "Queue"="Geometry"
            "RenderType"="Opaque"
        }
        LOD 100
        HLSLINCLUDE
        #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
        #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
        #include "ltc_input.hlsl"
        #include "forwardPass.hlsl"
        ENDHLSL
        Pass
        {
            Tags{"LightMode"="UniversalForward"}
            HLSLPROGRAM
            #pragma vertex ForwardPassVertex
            #pragma fragment ForwardPassFragment
            // make fog work
            #pragma multi_compile_fog
            
            ENDHLSL
        }
    }
}
