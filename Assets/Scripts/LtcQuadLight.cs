using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Unity.Mathematics;
using Unity.VisualScripting;
using UnityEngine;

public class LtcQuadLight : MonoBehaviour
{
    // Start is called before the first frame update
    public List<Texture2D> Textures;

    [HideInInspector]
    public Texture2DArray TextureArray;
    [HideInInspector]
    public ComputeBuffer VertPosBuffer;

    public static LtcQuadLight Instance;
    private Mesh m_mesh;
    
    void Awake()
    {
        Instance = this;
        m_mesh = GetComponent<MeshFilter>().sharedMesh;
        VertPosBuffer = new ComputeBuffer(4, Marshal.SizeOf(typeof(float3)));
        Refresh();
    }

    void Refresh()
    {
        TextureArray = new Texture2DArray(Textures[0].width,
            Textures[0].height, Textures.Count, Textures[0].format, true, false);
        // Apply settings
        TextureArray.filterMode = FilterMode.Bilinear;
        TextureArray.wrapMode = TextureWrapMode.Clamp;

        for (int i = 0; i < Textures.Count; i++)
        {
            Graphics.CopyTexture(Textures[i],0,TextureArray, i);
        }

        // Apply our changes
        TextureArray.Apply(false);
    }

    // Update is called once per frame
    void Update()
    {
        //VertPosBuffer.Dispose();
        var VertPos = new float3[4];
        /*var p0 = transform.localToWorldMatrix * new Vector4(0.5f, 0.5f, 0.0f, 1.0f);
        VertPos[0] = new float3(p0.x,p0.y,p0.z);
        
        var p1 = transform.localToWorldMatrix * new Vector4(-0.5f, 0.5f, 0.0f, 1.0f);
        VertPos[1] = new float3(p1.x,p1.y,p1.z);
        
        var p2 = transform.localToWorldMatrix * new Vector4(-0.5f, -0.5f, 0.0f, 1.0f);
        VertPos[2] = new float3(p2.x,p2.y,p2.z);
        
        var p3 = transform.localToWorldMatrix * new Vector4(-0.5f, -0.5f, 0.0f, 1.0f);
        VertPos[3] = new float3(p3.x,p3.y,p3.z);*/
        var p0 = transform.localToWorldMatrix * new Vector4(0.5f, 0.5f, 0.0f, 1.0f);
        VertPos[0] = new float3(p0.x,p0.y,p0.z);
        
        var p1 = transform.localToWorldMatrix * new Vector4(0.5f, -0.5f, 0.0f, 1.0f);
        VertPos[3] = new float3(p1.x,p1.y,p1.z);
        
        var p2 = transform.localToWorldMatrix * new Vector4(-0.5f, -0.5f, 0.0f, 1.0f);
        VertPos[2] = new float3(p2.x,p2.y,p2.z);
        
        var p3 = transform.localToWorldMatrix * new Vector4(-0.5f, 0.5f, 0.0f, 1.0f);
        VertPos[1] = new float3(p3.x,p3.y,p3.z);
        
        /*var p0 = transform.localToWorldMatrix * new Vector4(-0.5f, 0.5f, 0.0f, 1.0f);
        VertPos[0] = new float3(p0.x,p0.y,p0.z);
        
        var p1 = transform.localToWorldMatrix * new Vector4(-0.5f, -0.5f, 0.0f, 1.0f);
        VertPos[3] = new float3(p1.x,p1.y,p1.z);
        
        var p2 = transform.localToWorldMatrix * new Vector4(0.5f, -0.5f, 0.0f, 1.0f);
        VertPos[2] = new float3(p2.x,p2.y,p2.z);
        
        var p3 = transform.localToWorldMatrix * new Vector4(0.5f, 0.5f, 0.0f, 1.0f);
        VertPos[1] = new float3(p3.x,p3.y,p3.z);*/
        
        VertPosBuffer.SetData(VertPos);
    }
}
