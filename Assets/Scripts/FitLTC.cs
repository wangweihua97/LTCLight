/*﻿
using System;
using Unity.Mathematics;
using UnityEngine;

namespace DefaultNamespace
{
    public class FitLTC
    {
        // size of precomputed table (theta, alpha)
        static int N = 64;	//64
// number of samples used to compute the error during fitting
        static int Nsample = 32;	//32
// minimal roughness (avoid singularities)
        static float MIN_ALPHA = 0.00001f;

        static float pi = Mathf.PI;

// computes
// * the norm (albedo) of the BRDF
// * the average Schlick Fresnel value
// * the average direction of the BRDF
// * 其中norm、fresnel和averageDir都是乘了BRDF，并除了pdf的
void computeAvgTerms(BrdfGGX brdf,  float3 V,  float alpha,
    out float norm, out float fresnel,out float3 averageDir)
{
    norm = 0.0f;
    fresnel = 0.0f;
    averageDir = new float3(0, 0, 0);

    for (int j = 0; j < Nsample; ++j)	//在整个0~1的随机范围内去计算L，计算brdf的值，采样数越多，0~1范围内被覆盖的点越多，算出来的L以及brdf的值越精准（感觉随机生成U1、U2的值也可以）
    for (int i = 0; i < Nsample; ++i)
    {
         float U1 = (i + 0.5f)/Nsample;
         float U2 = (j + 0.5f)/Nsample;

        // sample
         float3 L = brdf.sample(V, alpha, U1, U2);

        // eval
        float pdf;
        float eval = brdf.eval(V, L, alpha,out pdf);

        if (pdf > 0)
        {
            float weight = eval / pdf;

            float3 H = math.normalize(V+L);

            // accumulate
            norm       += weight;	//norm存储的是论文Fresnel项附加材料中的nD
            fresnel    += weight * math.pow(1.0f - math.max(math.dot(V, H), 0.0f), 5.0f);		//fresnel存储的是论文Fresnel项附加材料中的fD
            averageDir += weight * L;
        }
    }

    norm    /= (float)(Nsample*Nsample);
    fresnel /= (float)(Nsample*Nsample);

    // clear y component, which should be zero with isotropic BRDFs
    averageDir.y = 0.0f;

    averageDir = math.normalize(averageDir);
}

// compute the error between the BRDF and the LTC
// using Multiple Importance Sampling
// 分别从BRDF和LTC里进行重要性采样，采出L，然后计算两种分布的值的差异，除以联合密度
static float computeError( LTC ltc,  BrdfGGX brdf,  float3 V,  float alpha)
{
    double error = 0.0;

    for (int j = 0; j < Nsample; ++j)
    for (int i = 0; i < Nsample; ++i)
    {
         float U1 = (i + 0.5f)/Nsample;
         float U2 = (j + 0.5f)/Nsample;

        // importance sample LTC
        {
            // sample
             float3 L = ltc.sample(U1, U2);

            float pdf_brdf;
            float eval_brdf = brdf.eval(V, L, alpha,out pdf_brdf);
            float eval_ltc = ltc.eval(L);
            float pdf_ltc = eval_ltc/ltc.magnitude;

            // error with MIS weight
            double error_ = MathF.Abs(eval_brdf - eval_ltc);
            error_ = error_*error_*error_;
            error += error_/(pdf_ltc + pdf_brdf);
        }

        // importance sample BRDF
        {
            // sample
             float3 L = brdf.sample(V, alpha, U1, U2);

            float pdf_brdf;
            float eval_brdf = brdf.eval(V, L, alpha,out pdf_brdf);
            float eval_ltc = ltc.eval(L);
            float pdf_ltc = eval_ltc/ltc.magnitude;

            // error with MIS weight
            double error_ = MathF.Abs(eval_brdf - eval_ltc);
            error_ = error_*error_*error_;
            error += error_/(pdf_ltc + pdf_brdf);
        }
    }

    return (float)error / (float)(Nsample*Nsample);
}

public class Fit_LTC
{
    BrdfGGX brdf;
    LTC ltc;
    bool isotropic;

    float3 V;
    float alpha;
    public Fit_LTC(LTC ltc_,  BrdfGGX brdf, bool isotropic_,  float3 V_, float alpha_)
         //brdf(brdf), V(V_), alpha(alpha_), isotropic(isotropic_)
    {
        ltc = ltc_;
        this.brdf = brdf;
        this.isotropic = isotropic_;
        V = V_;
        alpha = alpha_;
    }

    void update(float[] p)
    {
        float m11 = math.max(p[0], 1e-7f);
        float m22 = math.max(p[1], 1e-7f);
        float m13 = p[2];

        if (isotropic)
        {
            ltc.m11 = m11;
            ltc.m22 = m11;
            ltc.m13 = 0.0f;
        }
        else
        {
            ltc.m11 = m11;
            ltc.m22 = m22;
            ltc.m13 = m13;
        }
        ltc.update();
    }

    public float Test (float[] p)
    {
        update(p);
        return computeError(ltc, brdf, V, alpha);
    }

     
};
void mov(float3[] r,  float[] v, int dim)
{
    for (int i = 0; i < dim; ++i)
        r[i] = v[i];
}
float NelderMead(
    float[] pmin, float[] start, float delta, float tolerance, int maxIters, FUNC objectiveFn)
{
    // standard coefficients from Nelder-Mead
    const float reflect  = 1.0f;
    const float expand   = 2.0f;
    const float contract = 0.5f;
    const float shrink   = 0.5f;
    const int DIM = 3;

    const int NB_POINTS = DIM + 1;

    float3[] s = new float3[NB_POINTS];
    float[] f = new float[NB_POINTS];

    // initialise simplex
    mov(s, start, DIM);
    for (int i = 1; i < NB_POINTS; i++)
    {
        mov(s, start, DIM);
        s[i][i - 1] += delta;	//加0.05
    }

    // evaluate function at each point on simplex
    for (int i = 0; i < NB_POINTS; i++)
        f[i] = objectiveFn(s[i]);	//f[i]是LTC和BRDF在某个观察方向以及粗糙度下球面分布的差异，不同的f[i]对应不同的矩阵的3个元素m11、m22、m13

    int lo = 0, hi, nh;

    for (int j = 0; j < maxIters; j++)
    {
        // find lowest, highest and next highest，找出第几个f[i]是最小的、最大的和次大的（点NB_POINTS的序号）
        lo = hi = nh = 0;
        for (int i = 1; i < NB_POINTS; i++)
        {
            if (f[i] < f[lo])
                lo = i;
            if (f[i] > f[hi])
            {
                nh = hi;
                hi = i;
            }
            else if (f[i] > f[nh])
                nh = i;
        }

        // stop if we've reached the required tolerance level
        float a = fabsf(f[lo]);
        float b = fabsf(f[hi]);
        if (2.0f*fabsf(a - b) < (a + b)*tolerance)
            break;

        // compute centroid (excluding the worst point)
        point o;
        set(o, 0.0f, DIM);
        for (int i = 0; i < NB_POINTS; i++)
        {
            if (i == hi) continue;
            add(o, s[i], DIM);	//把除了最高误差之外的点对应的矩阵元素加到o里
        }

        for (int i = 0; i < DIM; i++)
            o[i] /= DIM;		//求均值

        // reflection
        float3 r;
        for (int i = 0; i < DIM; i++)
            r[i] = o[i] + reflect*(o[i] - s[hi][i]);	//点r中存储的是除了最高误差只玩的点的矩阵元素的均值，加上均值到最高误差的矩阵元素的差值，其实结果就是最高误差的矩阵元素

        float fr = objectiveFn(r);
        if (fr < f[nh])
        {
            if (fr < f[lo])
            {
                // expansion
                point e;
                for (int i = 0; i < DIM; i++)
                    e[i] = o[i] + expand*(o[i] - s[hi][i]);

                float fe = objectiveFn(e);
                if (fe < fr)
                {
                    mov(s[hi], e, DIM);
                    f[hi] = fe;
                    continue;
                }
            }

            mov(s[hi], r, DIM);
            f[hi] = fr;
            continue;
        }

        // contraction
        point c;
        for (int i = 0; i < DIM; i++)
            c[i] = o[i] - contract*(o[i] - s[hi][i]);

        float fc = objectiveFn(c);
        if (fc < f[hi])
        {
            mov(s[hi], c, DIM);
            f[hi] = fc;
            continue;
        }

        // reduction
        for (int k = 0; k < NB_POINTS; k++)
        {
            if (k == lo) continue;
            for (int i = 0; i < DIM; i++)
                s[k][i] = s[lo][i] + shrink*(s[k][i] - s[lo][i]);
            f[k] = objectiveFn(s[k]);
        }
    }

    // return best point and its value
    mov(pmin, s[lo], DIM);
    return f[lo];	//返回最低误差
}

// fit brute force
// refine first guess by exploring parameter space
void fit(LTC ltc,  BrdfGGX brdf,  float3 V,  float alpha,  float epsilon = 0.05f,  bool isotropic = false)
{
    float3 startFit = new float3( ltc.m11, ltc.m22, ltc.m13 );
    float3 resultFit;

    Fit_LTC fitter = new Fit_LTC(ltc, brdf, isotropic, V, alpha);

    // Find best-fit LTC lobe (scale, alphax, alphay)
    float error = NelderMead<3>(resultFit, startFit, epsilon, 1e-5f, 100, fitter);

    // Update LTC with best fitting values
    fitter.update(resultFit);
}

// fit data
// tab里存储的是LTC矩阵，tabMagFresnel里存储的是LTC的magnitude和fresnel
void fitTab(float3x3 tab, float2 tabMagFresnel,  int N,  BrdfGGX brdf)
{
    LTC ltc = new LTC();

    //  第一层循环枚举粗糙度，第二层循环枚举观察向量与平面的夹角从而得到与法线的夹角theta
    for (int a = N - 1; a >=     0; --a)	//粗糙度从1变到0
    for (int t =     0; t <= N - 1; ++t)	//theta从0变到1
    {
        // parameterised by sqrt(1 - cos(theta)), 即x = sqrt(1 - cos(theta))，所以cos(theta) = 1 - x^2
        float x = (float)t/(N - 1);
        float ct = 1.0f - x*x;
        float theta = math.min(1.57f, Mathf.Acos(ct));
        float3 V = new float3(Mathf.Sin(theta), 0, Mathf.Cos(theta));	//y分量是0，是因为在切线空间下算的，y轴垂直V和N所在的平面

        // alpha = roughness^2
        float roughness = a/(float)(N - 1);
        float alpha = math.max(roughness*roughness, MIN_ALPHA);	//alpha是平方粗糙度，也就是说在计算BRDF和LTC的值时用的是平方粗糙度


        float3 averageDir;
        computeAvgTerms(brdf, V, alpha, ltc.magnitude, ltc.fresnel, averageDir);

        bool isotropic;

        // 1.首先要猜测一个M矩阵
        // 初始化一个较为合适的半球体分布函数D0
        // 如果theta==0，则lobe是关于Z轴各向对称的
        if (t == 0)
        {
            ltc.X = float3(1, 0, 0);
            ltc.Y = float3(0, 1, 0);
            ltc.Z = float3(0, 0, 1);

            if (a == N - 1) // roughness = 1
            {
                ltc.m11 = 1.0f;
                ltc.m22 = 1.0f;
            }
            else // init with roughness of previous fit
            {
                ltc.m11 = tab[a + 1 + t*N][0][0];
                ltc.m22 = tab[a + 1 + t*N][1][1];
            }

            ltc.m13 = 0;
            ltc.update();	//初始化LTC矩阵

            isotropic = true;
        }
        // otherwise use previous configuration as first guess
        else
        {
            float3 L = averageDir;
            float3 T1(L.z, 0, -L.x);
            float3 T2(0, 1, 0);
            ltc.X = T1;
            ltc.Y = T2;
            ltc.Z = L;

            ltc.update();	//将LTC矩阵旋转到和向量L对齐

            isotropic = false;
        }

        // 2. fit (explore parameter space and refine first guess)通过单纯形法不断修正M矩阵使积分结果无限近似BRDF
        float epsilon = 0.05f;
        fit(ltc, brdf, V, alpha, epsilon, isotropic);

        // copy data
        tab[a + t*N] = ltc.M;	//按列存储，同一列上是theta在变
        tabMagFresnel[a + t*N][0] = ltc.magnitude;
        tabMagFresnel[a + t*N][1] = ltc.fresnel;

        // kill useless coefs in matrix，矩阵里只留下5个元素
        tab[a+t*N][0][1] = 0;
        tab[a+t*N][1][0] = 0;
        tab[a+t*N][2][1] = 0;
        tab[a+t*N][1][2] = 0;

        cout << tab[a+t*N][0][0] << "\t " << tab[a+t*N][1][0] << "\t " << tab[a+t*N][2][0] << endl;
        cout << tab[a+t*N][0][1] << "\t " << tab[a+t*N][1][1] << "\t " << tab[a+t*N][2][1] << endl;
        cout << tab[a+t*N][0][2] << "\t " << tab[a+t*N][1][2] << "\t " << tab[a+t*N][2][2] << endl;
        cout << endl;
    }
}

float sqr(float x)
{
    return x*x;
}

float G(float w, float s, float g)
{
    return -2.0f*sinf(w)*cosf(s)*cosf(g) + pi/2.0f - g + sinf(g)*cosf(g);
}

float H(float w, float s, float g)
{
    float sinsSq = sqr(sin(s));
    float cosgSq = sqr(cos(g));

    return cosf(w)*(cosf(g)*sqrtf(sinsSq - cosgSq) + sinsSq*asinf(cosf(g)/sinf(s)));
}

float ihemi(float w, float s)
{
    float g = asinf(cosf(s)/sinf(w));
    float sinsSq = sqr(sinf(s));

    if (w >= 0.0f && w <= (pi/2.0f - s))
        return pi*cosf(w)*sinsSq;

    if (w >= (pi/2.0f - s) && w < pi/2.0f)
        return pi*cosf(w)*sinsSq + G(w, s, g) - H(w, s, g);

    if (w >= pi/2.0f && w < (pi/2.0f + s))
        return G(w, s, g) + H(w, s, g);

    return 0.0f;
}

void genSphereTab(float* tabSphere, int N)
{
    for (int j = 0; j < N; ++j)
    for (int i = 0; i < N; ++i)
    {
         float U1 = float(i)/(N - 1);
         float U2 = float(j)/(N - 1);

        // z = cos(elevation angle)
        float z = 2.0f*U1 - 1.0f;

        // length of average dir., proportional to sin(sigma)^2
        float len = U2;

        float sigma = asinf(sqrtf(len));
        float omega = acosf(z);

        // compute projected (cosine-weighted) solid angle of spherical cap
        float value = 0.0f;

        if (sigma > 0.0f)
            value = ihemi(omega, sigma)/(pi*len);
        else
            value = std::max<float>(z, 0.0f);

        if (value != value)
            printf("nan!\n");

        tabSphere[i + j*N] = value;
    }
}

//把tab中存储的LTC矩阵的逆矩阵、tabMagFresnel中存储的LTC的Magnitude和Fresnel、以及tabSphere中存储的立体角存到纹理里
void packTab(
    vec4* tex1, vec4* tex2,
     mat3*  tab,
     vec2*  tabMagFresnel,
     float* tabSphere,
    int N)
{
    for (int i = 0; i < N*N; ++i)
    {
         mat3 m = tab[i];

        mat3 invM = inverse(m);

        // normalize by the middle element，除了矩阵的中间元素，可以少存一个元素（从5个变成4个）
        invM /= invM[1][1];

        // store the variable terms
        tex1[i].x = invM[0][0];
        tex1[i].y = invM[0][2];
        tex1[i].z = invM[2][0];
        tex1[i].w = invM[2][2];
        tex2[i].x = tabMagFresnel[i][0];
        tex2[i].y = tabMagFresnel[i][1];
        tex2[i].z = 0.0f; // unused
        tex2[i].w = tabSphere[i];
    }
}

int main(int argc, char* argv[])
{
    // BRDF to fit
    BrdfGGX brdf;
    //BrdfBeckmann brdf;
    //BrdfDisneyDiffuse brdf;

    // allocate data
    mat3*  tab = new mat3[N*N];
    vec2*  tabMagFresnel = new vec2[N*N];
    float* tabSphere = new float[N*N];

    // fit
    fitTab(tab, tabMagFresnel, N, brdf);

    // projected solid angle of a spherical cap, clipped to the horizon
    genSphereTab(tabSphere, N);

    // pack tables (texture representation)
    vec4* tex1 = new vec4[N*N];
    vec4* tex2 = new vec4[N*N];
    packTab(tex1, tex2, tab, tabMagFresnel, tabSphere, N);	//把tab中存储的LTC矩阵的逆矩阵、tabMagFresnel中存储的LTC的Magnitude和Fresnel、以及tabSphere中存储的立体角存到纹理里

    // export to C, MATLAB and DDS
    writeTabMatlab(tab, tabMagFresnel, N);
    writeTabC(tab, tabMagFresnel, N);
    writeDDS(tex1, tex2, N);
    writeJS(tex1, tex2, N);

    // spherical plots
    make_spherical_plots(brdf, tab, N);

    // delete data
    delete[] tab;
    delete[] tabMagFresnel;
    delete[] tabSphere;
    delete[] tex1;
    delete[] tex2;

    return 0;
}
    }
}*/