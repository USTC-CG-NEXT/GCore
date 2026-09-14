#include "GCore/Components/MaterialComponent.h"

#include <cstdio>
#include <filesystem>

RUZINO_NAMESPACE_OPEN_SCOPE
pxr::UsdShadeMaterial MaterialComponent::define_material(
    pxr::UsdStageRefPtr stage,
    pxr::SdfPath path) const
{
    auto material_path = get_material_path();

    auto material_shader_path =
        material_path.AppendPath(pxr::SdfPath("PBRShader"));
    pxr::UsdShadeMaterial material =
        pxr::UsdShadeMaterial::Define(stage, material_path);
    auto pbrShader = pxr::UsdShadeShader::Define(stage, material_shader_path);

    pbrShader.CreateIdAttr(pxr::VtValue(pxr::TfToken("UsdPreviewSurface")));
    material.CreateSurfaceOutput().ConnectToSource(
        pbrShader.ConnectableAPI(), pxr::TfToken("surface"));

    if (!material_object.textured || material_object.textures.empty()) {
        // Flat-color PBR: constant diffuse, no texture network, no primvar
        // plumbing needed.
        pbrShader
            .CreateInput(
                pxr::TfToken("diffuseColor"), pxr::SdfValueTypeNames->Color3f)
            .Set(
                pxr::GfVec3f(
                    material_object.base_color.r,
                    material_object.base_color.g,
                    material_object.base_color.b));
        pbrShader
            .CreateInput(
                pxr::TfToken("roughness"), pxr::SdfValueTypeNames->Float)
            .Set(material_object.roughness);
        return material;
    }

    auto material_stReader_path =
        material_path.AppendPath(pxr::SdfPath("stReader"));
    auto material_texture_path =
        material_path.AppendPath(pxr::SdfPath("diffuseTexture"));

    const std::string texture_name =
        std::string(material_object.textures[0].c_str());
    auto stReader = pxr::UsdShadeShader::Define(stage, material_stReader_path);
    stReader.CreateIdAttr(
        pxr::VtValue(pxr::TfToken("UsdPrimvarReader_float2")));

    auto diffuseTextureSampler =
        pxr::UsdShadeShader::Define(stage, material_texture_path);

    diffuseTextureSampler.CreateIdAttr(
        pxr::VtValue(pxr::TfToken("UsdUVTexture")));
    diffuseTextureSampler
        .CreateInput(pxr::TfToken("file"), pxr::SdfValueTypeNames->Asset)
        .Set(pxr::SdfAssetPath(texture_name));
    diffuseTextureSampler
        .CreateInput(pxr::TfToken("st"), pxr::SdfValueTypeNames->Float2)
        .ConnectToSource(stReader.ConnectableAPI(), pxr::TfToken("result"));
    diffuseTextureSampler.CreateOutput(
        pxr::TfToken("rgb"), pxr::SdfValueTypeNames->Float3);

    diffuseTextureSampler
        .CreateInput(pxr::TfToken("wrapS"), pxr::SdfValueTypeNames->Token)
        .Set(pxr::TfToken(material_object.wrap_mode.c_str()));

    diffuseTextureSampler
        .CreateInput(pxr::TfToken("wrapT"), pxr::SdfValueTypeNames->Token)
        .Set(pxr::TfToken(material_object.wrap_mode.c_str()));

    pbrShader
        .CreateInput(
            pxr::TfToken("diffuseColor"), pxr::SdfValueTypeNames->Color3f)
        .ConnectToSource(
            diffuseTextureSampler.ConnectableAPI(), pxr::TfToken("rgb"));

    pbrShader
        .CreateInput(pxr::TfToken("roughness"), pxr::SdfValueTypeNames->Float)
        .Set(material_object.roughness);

    if (material_object.alpha_cutout) {
        // Alpha-cutout: bind the same texture's alpha channel as opacity.
        // The path tracer's stochastic presence test (and the raster path)
        // drops samples below opacityThreshold, so transparent texels never
        // shade and shadow rays pass through — leaf-card / foliage style.
        auto material_opacity_texture_path =
            material_path.AppendPath(pxr::SdfPath("opacityTexture"));
        auto opacityTextureSampler =
            pxr::UsdShadeShader::Define(stage, material_opacity_texture_path);
        opacityTextureSampler.CreateIdAttr(
            pxr::VtValue(pxr::TfToken("UsdUVTexture")));
        opacityTextureSampler
            .CreateInput(pxr::TfToken("file"), pxr::SdfValueTypeNames->Asset)
            .Set(pxr::SdfAssetPath(texture_name));
        opacityTextureSampler
            .CreateInput(pxr::TfToken("st"), pxr::SdfValueTypeNames->Float2)
            .ConnectToSource(stReader.ConnectableAPI(), pxr::TfToken("result"));
        opacityTextureSampler.CreateOutput(
            pxr::TfToken("a"), pxr::SdfValueTypeNames->Float);
        opacityTextureSampler
            .CreateInput(pxr::TfToken("wrapS"), pxr::SdfValueTypeNames->Token)
            .Set(pxr::TfToken(material_object.wrap_mode.c_str()));
        opacityTextureSampler
            .CreateInput(pxr::TfToken("wrapT"), pxr::SdfValueTypeNames->Token)
            .Set(pxr::TfToken(material_object.wrap_mode.c_str()));

        pbrShader
            .CreateInput(pxr::TfToken("opacity"), pxr::SdfValueTypeNames->Float)
            .ConnectToSource(
                opacityTextureSampler.ConnectableAPI(), pxr::TfToken("a"));
        pbrShader
            .CreateInput(
                pxr::TfToken("opacityMode"), pxr::SdfValueTypeNames->Token)
            .Set(pxr::TfToken("mask"));
        pbrShader
            .CreateInput(
                pxr::TfToken("opacityThreshold"), pxr::SdfValueTypeNames->Float)
            .Set(material_object.opacity_threshold);
    }

    auto stInput = material.CreateInput(
        pxr::TfToken("frame:stPrimvarName"), pxr::SdfValueTypeNames->Token);
    stInput.Set(pxr::TfToken("UVMap"));

    stReader.CreateInput(pxr::TfToken("varname"), pxr::SdfValueTypeNames->Token)
        .ConnectToSource(stInput);

    return material;
}

pxr::SdfPath MaterialComponent::get_material_path() const
{
    if (!material_object.textured || material_object.textures.empty()) {
        // Color-derived name so equal colors share one material definition
        char buf[64];
        std::snprintf(
            buf,
            sizeof(buf),
            "flat_%d_%d_%d_Mat",
            int(material_object.base_color.r * 255.f + 0.5f),
            int(material_object.base_color.g * 255.f + 0.5f),
            int(material_object.base_color.b * 255.f + 0.5f));
        return pxr::SdfPath("/TexModel").AppendPath(pxr::SdfPath(buf));
    }
    auto texture_name = std::string(material_object.textures[0].c_str());
    std::filesystem::path p =
        std::filesystem::path(texture_name).replace_extension();
    auto file_name = "texture" + p.filename().string();
    auto material_path_root = pxr::SdfPath("/TexModel");
    auto material_path =
        material_path_root.AppendPath(pxr::SdfPath(file_name + "Mat"));
    return material_path;
}

RUZINO_NAMESPACE_CLOSE_SCOPE
