using System;
using UnityEngine;

namespace DefaultNamespace
{
    public class Ltc_Go : MonoBehaviour
    {
        private Material m_mat;
        
        private void Start()
        {
            m_mat = GetComponent<MeshRenderer>().sharedMaterial;
            //m_mat.SetTexture("_FilteredLightTexture" ,LtcQuadLight.Instance.TextureArray);
            m_mat.SetTexture("_FilteredLightTexture" ,Ltc_video.VideoList[0].TextureArray);
        }

        private void OnWillRenderObject()
        {
            //m_mat.SetTexture("_FilteredLightTexture" ,LtcQuadLight.Instance.TextureArray);
            //m_mat.SetBuffer("_PolygonalLightVertexPos",LtcQuadLight.Instance.VertPosBuffer);
            m_mat.SetTexture("_FilteredLightTexture" ,Ltc_video.VideoList[0].TextureArray);
            m_mat.SetBuffer("_PolygonalLightVertexPos",Ltc_video.VideoList[0].VertPosBuffer);
        }
    }
}