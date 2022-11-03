using System;
using UnityEngine;

namespace DefaultNamespace
{
    public class BrdfGGX
    {
        public  float eval(Vector3 V, Vector3 L,  float alpha,out float pdf)
        {
            if (V.z <= 0)
            {
                pdf = 0;
                return 0;
            }

            // masking
            float LambdaV = lambda(alpha, V.z);

            // shadowing
            float G2;
            if (L.z <= 0.0f)
                G2 = 0;
            else
            {
                float LambdaL = lambda(alpha, L.z);
                G2 = 1.0f/(1.0f + LambdaV + LambdaL);
            }

            // D
            Vector3 H = (V + L).normalized;
            float slopex = H.x/H.z;
            float slopey = H.y/H.z;
            float D = 1.0f / (1.0f + (slopex*slopex + slopey*slopey)/alpha/alpha);	//这个slopex*slopex + slopey*slopey其实等于(x^2+y^2)/z^2，结果就是tan(theta) * tan(theta)
            D = D*D;
            D = D/(3.14159f * alpha*alpha * H.z*H.z*H.z*H.z);

            pdf = MathF.Abs(D * H.z / 4.0f / Vector3.Dot(V, H));
            float res = D * G2 / 4.0f / V.z;		//这里实际上算的是BRDF * dot(N, L)，只是被BRDF分母里的dot(N, L)约掉了

            return res;
        }

        public  virtual Vector3 sample( Vector3 V, float alpha, float U1, float U2)
        {
            float phi = 2.0f*3.14159f * U1;
            float r = alpha*Mathf.Sqrt(U2/(1.0f - U2));
            Vector3 N = new Vector3(r*Mathf.Cos(phi), r*Mathf.Sin(phi), 1.0f).normalized;
            Vector3 L = -V + 2.0f * N * Vector3.Dot(N, V);
            return L;
        }
        private float lambda(float alpha, float cosTheta)
        {
            float a = 1.0f / alpha / Mathf.Tan(Mathf.Acos(cosTheta));
            return (cosTheta < 1.0f) ? 0.5f * (-1.0f + Mathf.Sqrt(1.0f + 1.0f/a/a)) : 0.0f;    
        }
    }
}