
cbuffer SearchParams : register(b0)
{
    float radius;
};

RaytracingAccelerationStructure SceneBVH : register(t0);
StructuredBuffer<float3> positions : register(t1);

RWStructuredBuffer<uint2> Pairs : register(u0);
RWStructuredBuffer<uint> PairsCount : register(u1);

struct HitInfo {
    float4 ShadedColorAndHitT : SHADED_COLOR_AND_HIT_T;
};

struct Attributes {
    float2 uv;
};

[shader("raygeneration")]
void RayGen()
{
    //// Trace the ray
    HitInfo payload;
    RayDesc ray;

    ray.TMin = radius / 4;
    ray.TMax = radius / 2;
    ray.Direction = float3(0, 0, 1);
    ray.Origin = positions[DispatchRaysIndex().x];

    TraceRay(SceneBVH, RAY_FLAG_NONE, 0xFF, 0, 0, 0, ray, payload);
}

[shader("anyhit")]
void AnyHit(
    inout HitInfo payload: SV_RayPayload,
    Attributes attrib: SV_IntersectionAttributes)
{
    uint id;
    if (DispatchRaysIndex().x != PrimitiveIndex()) {
        InterlockedAdd(PairsCount[0], 1u, id);
        Pairs[id] = int2(DispatchRaysIndex().x, PrimitiveIndex());
    }
    IgnoreHit();
}

// ---[ Closest Hit Shader ]---

[shader("closesthit")]
void ClosestHit(
    inout HitInfo payload: SV_RayPayload,
    Attributes attrib: SV_IntersectionAttributes)
{
}

[shader("intersection")]
void Intersection()
{
    Attributes attribute;
    attribute.uv = float2(PrimitiveIndex() / 100.f, 0);

    // PairsCount[0] = 114514;

    ReportHit(RayTCurrent(), 0, attribute);
}

// ---[ Miss Shader ]---

[shader("miss")]
void Miss(inout HitInfo payload: SV_RayPayload)
{
}
