using System;
using Unity.Mathematics;
using UnityEngine;

namespace DefaultNamespace
{
    public class LTC
    {
        // lobe magnitude
        public float magnitude;

        // Average Schlick Fresnel term
        public float fresnel;

        // parametric representation
        public float m11, m22, m13;
        public float3 X, Y, Z;

        // matrix representation
        public float3x3 M;
        public float3x3 invM;
        public float detM;

        public LTC()
        {
            magnitude = 1;
            fresnel = 1;
            m11 = 1;
            m22 = 1;
            m13 = 0;
            X = new float3(1, 0, 0);
            Y = new float3(0, 1, 0);
            Z = new float3(0, 0, 1);
            update();
        }

        public void update() // compute matrix from parameters
        {
            M = math.mul(new float3x3(X, Y, Z) ,
                         new float3x3(m11, 0, 0,
                             0, m22, 0,
                             m13, 0, 1));
            invM = math.inverse(M);
            detM = MathF.Abs(math.determinant(M));
        }

        public float eval( float3 L)
        {
            float3 Loriginal = math.normalize(math.mul(invM , L));
            
            float3 L_ =math.mul( M , Loriginal);

            float l = math.length(L_);
            float Jacobian = detM / (l*l*l);

            float D = 1.0f / 3.14159f * math.max(0.0f, Loriginal.z); 
		
            float res = magnitude * D / Jacobian;
            return res;
        }

        public float3 sample(float U1, float U2)
        {
             float theta = Mathf.Acos(Mathf.Sqrt(U1));
             float phi = 2.0f*3.14159f * U2;
             float3 L = math.normalize(math.mul(M , new float3(Mathf.Sin(theta)*Mathf.Cos(phi), Mathf.Sin(theta)*Mathf.Sin(phi), Mathf.Cos(theta))));	//把余弦重要性采样得到的向量乘上矩阵M，变成BRDF重要性采样的向量
            return L;
        }
    }
}