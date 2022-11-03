struct SRay
{
    float3 m_Origin;
    float3 m_Dir;
};


void clipPolygonalLightByZPlane(inout float3 vioPolygonalLightVertexPos[5], out int voPolygonalLightVertexNumAfterClipping)
{
	int Flag4HowManyVertexAboveZPlane = 0;
	if(vioPolygonalLightVertexPos[0].z > 0.0) Flag4HowManyVertexAboveZPlane += 1;
	if(vioPolygonalLightVertexPos[1].z > 0.0) Flag4HowManyVertexAboveZPlane += 2;
	if(vioPolygonalLightVertexPos[2].z > 0.0) Flag4HowManyVertexAboveZPlane += 4;
	if(vioPolygonalLightVertexPos[3].z > 0.0) Flag4HowManyVertexAboveZPlane += 8;

	voPolygonalLightVertexNumAfterClipping = 0;
	if(Flag4HowManyVertexAboveZPlane == 0)
	{
	}
	else if(Flag4HowManyVertexAboveZPlane == 1)
	{
		voPolygonalLightVertexNumAfterClipping = 3;
		vioPolygonalLightVertexPos[1] = -vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[0] + vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[1];
		vioPolygonalLightVertexPos[2] = -vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[0] + vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[3];
	}
	else if(Flag4HowManyVertexAboveZPlane == 2)
	{
		voPolygonalLightVertexNumAfterClipping = 3;
		vioPolygonalLightVertexPos[0] = -vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[1] + vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[0];
		vioPolygonalLightVertexPos[2] = -vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[1] + vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[2];
	}
	else if(Flag4HowManyVertexAboveZPlane == 3)
	{
		voPolygonalLightVertexNumAfterClipping = 4;
		vioPolygonalLightVertexPos[2] = -vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[1] + vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[2];
		vioPolygonalLightVertexPos[3] = -vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[0] + vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[3];
	}
	else if(Flag4HowManyVertexAboveZPlane == 4)
	{
		voPolygonalLightVertexNumAfterClipping = 3;
		vioPolygonalLightVertexPos[0] = -vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[2] + vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[3];
		vioPolygonalLightVertexPos[1] = -vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[2] + vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[1];
	}
	else if(Flag4HowManyVertexAboveZPlane == 5)
	{
		voPolygonalLightVertexNumAfterClipping = 0;
	}
	else if(Flag4HowManyVertexAboveZPlane == 6)
	{
		voPolygonalLightVertexNumAfterClipping = 4;
		vioPolygonalLightVertexPos[0] = -vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[1] + vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[0];
		vioPolygonalLightVertexPos[3] = -vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[2] + vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[3];
	}
	else if(Flag4HowManyVertexAboveZPlane == 7)
	{
		voPolygonalLightVertexNumAfterClipping = 5;
		vioPolygonalLightVertexPos[4] = -vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[0] + vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[3];
		vioPolygonalLightVertexPos[3] = -vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[2] + vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[3];
	}
	else if(Flag4HowManyVertexAboveZPlane == 8)
	{
		voPolygonalLightVertexNumAfterClipping = 3;
		vioPolygonalLightVertexPos[0] = -vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[3] + vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[0];
		vioPolygonalLightVertexPos[1] = -vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[3] + vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[2];
		vioPolygonalLightVertexPos[2] = vioPolygonalLightVertexPos[3];
	}
	else if(Flag4HowManyVertexAboveZPlane == 9)
	{
		voPolygonalLightVertexNumAfterClipping = 4;
		vioPolygonalLightVertexPos[1] = -vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[0] + vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[1];
		vioPolygonalLightVertexPos[2] = -vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[3] + vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[2];
	}
	else if(Flag4HowManyVertexAboveZPlane == 10)
	{
		voPolygonalLightVertexNumAfterClipping = 0;
	}
	else if(Flag4HowManyVertexAboveZPlane == 11)
	{
		voPolygonalLightVertexNumAfterClipping = 5;
		vioPolygonalLightVertexPos[4] = vioPolygonalLightVertexPos[3];
		vioPolygonalLightVertexPos[3] = -vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[3] + vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[2];
		vioPolygonalLightVertexPos[2] = -vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[1] + vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[2];
	}
	else if(Flag4HowManyVertexAboveZPlane == 12)
	{
		voPolygonalLightVertexNumAfterClipping = 4;
		vioPolygonalLightVertexPos[1] = -vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[2] + vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[1];
		vioPolygonalLightVertexPos[0] = -vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[3] + vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[0];
	}
	else if(Flag4HowManyVertexAboveZPlane == 13)
	{
		voPolygonalLightVertexNumAfterClipping = 5;
		vioPolygonalLightVertexPos[4] = vioPolygonalLightVertexPos[3];
		vioPolygonalLightVertexPos[3] = vioPolygonalLightVertexPos[2];
		vioPolygonalLightVertexPos[2] = -vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[2] + vioPolygonalLightVertexPos[2].z * vioPolygonalLightVertexPos[1];
		vioPolygonalLightVertexPos[1] = -vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[0] + vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[1];
	}
	else if(Flag4HowManyVertexAboveZPlane == 14)
	{
		voPolygonalLightVertexNumAfterClipping = 5;
		vioPolygonalLightVertexPos[4] = -vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[3] + vioPolygonalLightVertexPos[3].z * vioPolygonalLightVertexPos[0];
		vioPolygonalLightVertexPos[0] = -vioPolygonalLightVertexPos[0].z * vioPolygonalLightVertexPos[1] + vioPolygonalLightVertexPos[1].z * vioPolygonalLightVertexPos[0];
	}
	else if(Flag4HowManyVertexAboveZPlane == 15)
	{
		voPolygonalLightVertexNumAfterClipping = 4;
	}

	if(voPolygonalLightVertexNumAfterClipping == 3)
		vioPolygonalLightVertexPos[3] = vioPolygonalLightVertexPos[0];
	if(voPolygonalLightVertexNumAfterClipping == 4)
		vioPolygonalLightVertexPos[4] = vioPolygonalLightVertexPos[0];
}

float3 integrateEdge(float3 vVertex1, float3 vVertex2)
{
	float CosTheta = dot(vVertex1, vVertex2);
	float Theta = acos(CosTheta);
	return cross(vVertex1, vVertex2) * ((Theta > 0.001) ? Theta / sin(Theta) : 1.0);
}

bool rayPlaneIntersect(SRay vRay, float4 vPlane, out float voDistanceFromRayOrigin2Plane)
{
    voDistanceFromRayOrigin2Plane = -dot(vPlane, float4(vRay.m_Origin, 1.0)) / dot(vPlane.xyz, vRay.m_Dir);
    return voDistanceFromRayOrigin2Plane > 0.0;
}

float3 fecthFilteredLightTexture(float3 vPolygonalLightVertexPos[4], float3 vLooupVector)
{
	float3 V1 = vPolygonalLightVertexPos[1] - vPolygonalLightVertexPos[0];
	float3 V2 = vPolygonalLightVertexPos[3] - vPolygonalLightVertexPos[0];
	float3 PlaneOrtho = cross(V1, V2);
	float PlaneAreaSquared = dot(PlaneOrtho, PlaneOrtho);
	
	SRay Ray;
	Ray.m_Origin = float3(0,0,0);
	Ray.m_Dir = vLooupVector;
	float4 Plane = float4(PlaneOrtho, -dot(PlaneOrtho, vPolygonalLightVertexPos[0]));
	float Distance2Plane;
	rayPlaneIntersect(Ray, Plane, Distance2Plane);
	
	float3 P = Distance2Plane * Ray.m_Dir - vPolygonalLightVertexPos[0];

	float Dot_V1_V2 = dot(V1, V2);
	float Inv_Dot_V1_V1 = 1.0 / dot(V1, V1);
	float3  V2_ = V2 - V1 * Dot_V1_V2 * Inv_Dot_V1_V1;
	float2  UV;
	UV.y = dot(V2_, P) / dot(V2_, V2_);
	UV.x = dot(V1, P) * Inv_Dot_V1_V1 - Dot_V1_V2 * Inv_Dot_V1_V1 * UV.y;
	UV = float2(1 - UV.y, 1 - UV.x);
	
	//test
	UV = float2(UV.y, UV.x);
	
	float Distance = abs(Distance2Plane) / pow(PlaneAreaSquared, 0.25);
	
	float Lod = log(2048.0 * Distance) / log(3.0);

	float LodA = floor(Lod);
	float LodB = ceil(Lod);
	float t = Lod - LodA;
	float3 ColorA = SAMPLE_TEXTURE2D_ARRAY(_FilteredLightTexture ,sampler_BaseMap, UV, LodA).rgb;
	float3 ColorB = SAMPLE_TEXTURE2D_ARRAY(_FilteredLightTexture ,sampler_BaseMap, UV, LodB).rgb;
	return lerp(ColorA, ColorB, t);
}


float3 integrateLTC(float3 vNormal, float3 vViewDir, float3 vFragPos, float3x3 vLTCMatrix, StructuredBuffer<float3> vPolygonalLightVertexPos, bool vTwoSided)
{	
	//着色点上的切线空间正交基
	float3 Tangent = normalize(vViewDir - (vNormal * dot(vViewDir, vNormal)));
	float3 Bitangent = cross(vNormal, Tangent);
	//float3 Bitangent = normalize(cross(vNormal, Tangent));

	//将变换矩阵转换到切线空间
	//vLTCMatrix = mul(vLTCMatrix ,transpose(float3x3(Tangent, Bitangent, vNormal)));
	vLTCMatrix = mul(vLTCMatrix ,float3x3(Tangent, Bitangent, vNormal));
	//vLTCMatrix = mul(float3x3(Tangent, Bitangent, vNormal) ,vLTCMatrix);
	//多边形顶点(因为被平面z=0裁剪以后，四边形可能变成5个顶点)
	float3 PolygonalLightVertexPosInTangentSpace[5], PolygonalLightVertexPosBeforeClipping[4];
	PolygonalLightVertexPosInTangentSpace[0] = mul(vLTCMatrix , (vPolygonalLightVertexPos[0] - vFragPos));
	PolygonalLightVertexPosBeforeClipping[0] = PolygonalLightVertexPosInTangentSpace[0];
	PolygonalLightVertexPosInTangentSpace[1] = mul(vLTCMatrix , (vPolygonalLightVertexPos[1] - vFragPos));
	PolygonalLightVertexPosBeforeClipping[1] = PolygonalLightVertexPosInTangentSpace[1];
	PolygonalLightVertexPosInTangentSpace[2] = mul(vLTCMatrix , (vPolygonalLightVertexPos[2] - vFragPos));
	PolygonalLightVertexPosBeforeClipping[2] = PolygonalLightVertexPosInTangentSpace[2];
	PolygonalLightVertexPosInTangentSpace[3] = mul(vLTCMatrix , (vPolygonalLightVertexPos[3] - vFragPos));
	
	/*PolygonalLightVertexPosInTangentSpace[0] = mul((vPolygonalLightVertexPos[0] - vFragPos) ,vLTCMatrix);
	PolygonalLightVertexPosBeforeClipping[0] = PolygonalLightVertexPosInTangentSpace[0];
	PolygonalLightVertexPosInTangentSpace[1] = mul((vPolygonalLightVertexPos[1] - vFragPos) ,vLTCMatrix);
	PolygonalLightVertexPosBeforeClipping[1] = PolygonalLightVertexPosInTangentSpace[1];
	PolygonalLightVertexPosInTangentSpace[2] = mul((vPolygonalLightVertexPos[2] - vFragPos) ,vLTCMatrix);
	PolygonalLightVertexPosBeforeClipping[2] = PolygonalLightVertexPosInTangentSpace[2];
	PolygonalLightVertexPosInTangentSpace[3] = mul((vPolygonalLightVertexPos[3] - vFragPos) ,vLTCMatrix);*/
	
	PolygonalLightVertexPosBeforeClipping[3] = PolygonalLightVertexPosInTangentSpace[3];
	PolygonalLightVertexPosInTangentSpace[4] = PolygonalLightVertexPosInTangentSpace[3];

	int PolygonalLightVertexNumAfterClipping;
	clipPolygonalLightByZPlane(PolygonalLightVertexPosInTangentSpace, PolygonalLightVertexNumAfterClipping);
	
	if(PolygonalLightVertexNumAfterClipping == 0)
		return float3(0.0,0.0,0.0);
	
	//把裁剪后的多边形投影到球面上（也就是对每个顶点坐标向量归一化）
	PolygonalLightVertexPosInTangentSpace[0] = normalize(PolygonalLightVertexPosInTangentSpace[0]);
	PolygonalLightVertexPosInTangentSpace[1] = normalize(PolygonalLightVertexPosInTangentSpace[1]);
	PolygonalLightVertexPosInTangentSpace[2] = normalize(PolygonalLightVertexPosInTangentSpace[2]);
	PolygonalLightVertexPosInTangentSpace[3] = normalize(PolygonalLightVertexPosInTangentSpace[3]);
	PolygonalLightVertexPosInTangentSpace[4] = normalize(PolygonalLightVertexPosInTangentSpace[4]);

	//用累加边的公式来求解积分
	float3 VSum = float3(0.0,0.0,0.0);
	VSum += integrateEdge(PolygonalLightVertexPosInTangentSpace[0], PolygonalLightVertexPosInTangentSpace[1]);
	VSum += integrateEdge(PolygonalLightVertexPosInTangentSpace[1], PolygonalLightVertexPosInTangentSpace[2]);
	VSum += integrateEdge(PolygonalLightVertexPosInTangentSpace[2], PolygonalLightVertexPosInTangentSpace[3]);
	if(PolygonalLightVertexNumAfterClipping >= 4)
		VSum += integrateEdge(PolygonalLightVertexPosInTangentSpace[3], PolygonalLightVertexPosInTangentSpace[4]);
	if(PolygonalLightVertexNumAfterClipping == 5)
		VSum += integrateEdge(PolygonalLightVertexPosInTangentSpace[4], PolygonalLightVertexPosInTangentSpace[0]);

	float Sum = vTwoSided > 0.1 ? abs(VSum.z) : max(VSum.z, 0.0);

	float3 LookupVector = normalize(VSum);
	float3 TextureLight = fecthFilteredLightTexture(PolygonalLightVertexPosBeforeClipping, LookupVector);

	return float3(Sum,Sum,Sum) * TextureLight;
	//return float3(VSum.z,VSum.z,VSum.z);
	//return float3(Sum,Sum,Sum);
}

float2 LTC_Coords(float vCosTheta, float vRoughness)
{
	float Theta = acos(vCosTheta);
    float2 Coords = float2(vRoughness, Theta/(0.5 * PI));
    const float LUT_SIZE = 32.0;
    // scale and bias coordinates, for correct filtered lookup
    Coords = Coords * (LUT_SIZE - 1.0)/LUT_SIZE + 0.5/LUT_SIZE;
    return Coords;
}