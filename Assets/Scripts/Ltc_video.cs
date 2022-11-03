using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Experimental.Rendering;
using UnityEngine.Rendering;
using UnityEngine.Video;

namespace DefaultNamespace
{
    public class Ltc_video : MonoBehaviour
    {
        public static List<Ltc_video> VideoList;
        public Material BlurMat;
        public float BlurSpread = 1.0f;
        [HideInInspector]
        public RenderTexture TextureArray;
        [HideInInspector]
        public ComputeBuffer VertPosBuffer;

        private VideoPlayer m_videoPlayer;

        private Material m_mat;

        private RenderTexture m_rt;
        //private List<RenderTexture> m_rts;

        private void Awake()
        {
            if(VideoList == null)
                VideoList = new List<Ltc_video>();
            VideoList.Add(this);
            m_videoPlayer = GetComponent<VideoPlayer>();
            VertPosBuffer = new ComputeBuffer(4, Marshal.SizeOf(typeof(float3)));
            m_mat = GetComponent<MeshRenderer>().sharedMaterial;
            m_rt = new RenderTexture(256,256,0,GraphicsFormat.R8G8B8A8_UNorm);
            m_mat.SetTexture("_BaseMap" ,m_rt);
            //m_rts = new List<RenderTexture>();
        }

        private void Start()
        {
            RenderTexture source = m_videoPlayer.targetTexture;
            RenderTextureDescriptor d = new RenderTextureDescriptor(source.width, source.height, RenderTextureFormat.ARGB32);
            d.dimension = TextureDimension.Tex2DArray;
            d.volumeDepth = 6; // We will have 2 slices (I.E. 2 textures) in our texture array.

            TextureArray = new RenderTexture(d);
            TextureArray.Create();
            TextureArray.name = "LTC Texture Array";
            //RenderTexture s = new RenderTexture();
            //Graphics.CopyTexture(s,0,TextureArray, 0);
        }

        private void Update()
        {
            var VertPos = new float3[4];
            var p0 = transform.localToWorldMatrix * new Vector4(0.5f, 0.5f, 0.0f, 1.0f);
            VertPos[0] = new float3(p0.x,p0.y,p0.z);
        
            var p1 = transform.localToWorldMatrix * new Vector4(0.5f, -0.5f, 0.0f, 1.0f);
            VertPos[3] = new float3(p1.x,p1.y,p1.z);
        
            var p2 = transform.localToWorldMatrix * new Vector4(-0.5f, -0.5f, 0.0f, 1.0f);
            VertPos[2] = new float3(p2.x,p2.y,p2.z);
        
            var p3 = transform.localToWorldMatrix * new Vector4(-0.5f, 0.5f, 0.0f, 1.0f);
            VertPos[1] = new float3(p3.x,p3.y,p3.z);
        
            VertPosBuffer.SetData(VertPos);
            RenderTexture source = m_videoPlayer.targetTexture;
            List<RenderTexture> rts = new List<RenderTexture>();
            Graphics.CopyTexture(source, 0,
                TextureArray, 0);
            for (int index = 0; index < 5; index++)
            {
                RenderTexture cur = source;
                for (int i = 0; i < 4; i++) {
                    BlurMat.SetFloat("_BlurSize", 1.0f + i * BlurSpread);
 
                    RenderTexture buffer0 = RenderTexture.GetTemporary(cur.width, cur.height, 0 ,source.format);
                    RenderTexture buffer1 = RenderTexture.GetTemporary(cur.width, cur.height, 0 ,source.format);
                    // Render the vertical pass
                    Graphics.Blit(cur, buffer0, BlurMat, 0);
                    // Render the horizontal pass
                    Graphics.Blit(buffer0, buffer1, BlurMat, 1);
 
                    // 空出buffer1为下一次循环准备
                    RenderTexture.ReleaseTemporary(buffer0);
                    if(i != 0)
                        RenderTexture.ReleaseTemporary(cur);
                    cur = buffer1;
                }
                rts.Add(cur);
                Graphics.CopyTexture(cur, 0,
                    TextureArray, index + 1);
                source = cur;
                
            }
            foreach (var rt in rts)
            {
                RenderTexture.ReleaseTemporary(rt);
            }
            m_mat.SetTexture("_BaseMap" ,m_videoPlayer.targetTexture);
            // Apply our changes
            
        }
    }
}