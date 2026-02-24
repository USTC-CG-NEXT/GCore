#pragma once

#include <cmath>
#include <cstdint>
#include <glm/glm.hpp>

RUZINO_NAMESPACE_OPEN_SCOPE

namespace Noise {

inline uint32_t hash(uint32_t x)
{
    x ^= x >> 17;
    x *= 0xed5ad4bb;
    x ^= x >> 11;
    x *= 0xac4c1b51;
    x ^= x >> 15;
    x *= 0x31848bab;
    x ^= x >> 14;
    return x;
}

inline uint32_t hash(uint32_t x, uint32_t y)
{
    return hash(x ^ hash(y));
}

inline uint32_t hash(uint32_t x, uint32_t y, uint32_t z)
{
    return hash(x ^ hash(y) ^ hash(z));
}

inline float hash_to_float(uint32_t h)
{
    return static_cast<float>(h) / static_cast<float>(0xffffffff);
}

inline float random_float(uint32_t seed)
{
    return hash_to_float(hash(seed));
}

inline float random_float(uint32_t x, uint32_t y)
{
    return hash_to_float(hash(x, y));
}

inline float random_float(uint32_t x, uint32_t y, uint32_t z)
{
    return hash_to_float(hash(x, y, z));
}

inline float random_float(const glm::vec2& p, uint32_t seed = 0)
{
    return hash_to_float(hash(
        static_cast<uint32_t>(std::floor(p.x * 8192.0f)) + seed,
        static_cast<uint32_t>(std::floor(p.y * 8192.0f))));
}

inline float random_float(const glm::vec3& p, uint32_t seed = 0)
{
    return hash_to_float(hash(
        static_cast<uint32_t>(std::floor(p.x * 8192.0f)) + seed,
        static_cast<uint32_t>(std::floor(p.y * 8192.0f)),
        static_cast<uint32_t>(std::floor(p.z * 8192.0f))));
}

inline glm::vec2 random_vec2(const glm::vec2& p, uint32_t seed = 0)
{
    float x = random_float(p, seed);
    float y = random_float(p, seed + 1);
    return glm::vec2(x, y);
}

inline glm::vec3 random_vec3(const glm::vec3& p, uint32_t seed = 0)
{
    float x = random_float(p, seed);
    float y = random_float(p, seed + 1);
    float z = random_float(p, seed + 2);
    return glm::vec3(x, y, z);
}

inline float value_noise_2d(const glm::vec2& p, uint32_t seed = 0)
{
    glm::vec2 i = glm::floor(p);
    glm::vec2 f = glm::fract(p);

    float a = random_float(
        static_cast<uint32_t>(i.x), static_cast<uint32_t>(i.y) + seed);
    float b = random_float(
        static_cast<uint32_t>(i.x + 1), static_cast<uint32_t>(i.y) + seed);
    float c = random_float(
        static_cast<uint32_t>(i.x), static_cast<uint32_t>(i.y + 1) + seed);
    float d = random_float(
        static_cast<uint32_t>(i.x + 1), static_cast<uint32_t>(i.y + 1) + seed);

    glm::vec2 u = f * f * (3.0f - 2.0f * f);

    return glm::mix(glm::mix(a, b, u.x), glm::mix(c, d, u.x), u.y);
}

inline float value_noise_3d(const glm::vec3& p, uint32_t seed = 0)
{
    glm::vec3 i = glm::floor(p);
    glm::vec3 f = glm::fract(p);

    float a = random_float(
        static_cast<uint32_t>(i.x),
        static_cast<uint32_t>(i.y),
        static_cast<uint32_t>(i.z) + seed);
    float b = random_float(
        static_cast<uint32_t>(i.x + 1),
        static_cast<uint32_t>(i.y),
        static_cast<uint32_t>(i.z) + seed);
    float c = random_float(
        static_cast<uint32_t>(i.x),
        static_cast<uint32_t>(i.y + 1),
        static_cast<uint32_t>(i.z) + seed);
    float d = random_float(
        static_cast<uint32_t>(i.x + 1),
        static_cast<uint32_t>(i.y + 1),
        static_cast<uint32_t>(i.z) + seed);
    float e = random_float(
        static_cast<uint32_t>(i.x),
        static_cast<uint32_t>(i.y),
        static_cast<uint32_t>(i.z + 1) + seed);
    float f1 = random_float(
        static_cast<uint32_t>(i.x + 1),
        static_cast<uint32_t>(i.y),
        static_cast<uint32_t>(i.z + 1) + seed);
    float g = random_float(
        static_cast<uint32_t>(i.x),
        static_cast<uint32_t>(i.y + 1),
        static_cast<uint32_t>(i.z + 1) + seed);
    float h = random_float(
        static_cast<uint32_t>(i.x + 1),
        static_cast<uint32_t>(i.y + 1),
        static_cast<uint32_t>(i.z + 1) + seed);

    glm::vec3 u = f * f * (3.0f - 2.0f * f);

    float x1 = glm::mix(a, b, u.x);
    float x2 = glm::mix(c, d, u.x);
    float x3 = glm::mix(e, f1, u.x);
    float x4 = glm::mix(g, h, u.x);

    float y1 = glm::mix(x1, x2, u.y);
    float y2 = glm::mix(x3, x4, u.y);

    return glm::mix(y1, y2, u.z);
}

inline glm::vec2 grad2d(uint32_t hash_val)
{
    static const glm::vec2 gradients[] = {
        glm::vec2(1.0f, 0.0f),        glm::vec2(-1.0f, 0.0f),
        glm::vec2(0.0f, 1.0f),        glm::vec2(0.0f, -1.0f),
        glm::vec2(0.7071f, 0.7071f),  glm::vec2(-0.7071f, 0.7071f),
        glm::vec2(0.7071f, -0.7071f), glm::vec2(-0.7071f, -0.7071f)
    };
    return gradients[hash_val & 7];
}

inline glm::vec3 grad3d(uint32_t hash_val)
{
    static const glm::vec3 gradients[] = {
        glm::vec3(1.0f, 1.0f, 0.0f),  glm::vec3(-1.0f, 1.0f, 0.0f),
        glm::vec3(1.0f, -1.0f, 0.0f), glm::vec3(-1.0f, -1.0f, 0.0f),
        glm::vec3(1.0f, 0.0f, 1.0f),  glm::vec3(-1.0f, 0.0f, 1.0f),
        glm::vec3(1.0f, 0.0f, -1.0f), glm::vec3(-1.0f, 0.0f, -1.0f),
        glm::vec3(0.0f, 1.0f, 1.0f),  glm::vec3(0.0f, -1.0f, 1.0f),
        glm::vec3(0.0f, 1.0f, -1.0f), glm::vec3(0.0f, -1.0f, -1.0f),
        glm::vec3(1.0f, 1.0f, 0.0f),  glm::vec3(-1.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, -1.0f, 1.0f), glm::vec3(0.0f, -1.0f, -1.0f)
    };
    return gradients[hash_val & 15];
}

inline float perlin_noise_2d(const glm::vec2& p, uint32_t seed = 0)
{
    glm::vec2 i = glm::floor(p);
    glm::vec2 f = glm::fract(p);

    glm::vec2 u = f * f * (3.0f - 2.0f * f);

    float a = glm::dot(
        grad2d(hash(
            static_cast<uint32_t>(i.x), static_cast<uint32_t>(i.y) + seed)),
        f);
    float b = glm::dot(
        grad2d(hash(
            static_cast<uint32_t>(i.x + 1), static_cast<uint32_t>(i.y) + seed)),
        f - glm::vec2(1.0f, 0.0f));
    float c = glm::dot(
        grad2d(hash(
            static_cast<uint32_t>(i.x), static_cast<uint32_t>(i.y + 1) + seed)),
        f - glm::vec2(0.0f, 1.0f));
    float d = glm::dot(
        grad2d(hash(
            static_cast<uint32_t>(i.x + 1),
            static_cast<uint32_t>(i.y + 1) + seed)),
        f - glm::vec2(1.0f, 1.0f));

    return glm::mix(glm::mix(a, b, u.x), glm::mix(c, d, u.x), u.y) * 0.5f +
           0.5f;
}

inline float perlin_noise_3d(const glm::vec3& p, uint32_t seed = 0)
{
    glm::vec3 i = glm::floor(p);
    glm::vec3 f = glm::fract(p);

    glm::vec3 u = f * f * (3.0f - 2.0f * f);

    uint32_t ix = static_cast<uint32_t>(i.x);
    uint32_t iy = static_cast<uint32_t>(i.y);
    uint32_t iz = static_cast<uint32_t>(i.z);

    auto dot_grad = [&](int dx, int dy, int dz, const glm::vec3& frac) {
        glm::vec3 g = grad3d(hash(ix + dx, iy + dy, iz + dz + seed));
        return glm::dot(g, frac - glm::vec3(dx, dy, dz));
    };

    float x1 = glm::mix(dot_grad(0, 0, 0, f), dot_grad(1, 0, 0, f), u.x);
    float x2 = glm::mix(dot_grad(0, 1, 0, f), dot_grad(1, 1, 0, f), u.x);
    float x3 = glm::mix(dot_grad(0, 0, 1, f), dot_grad(1, 0, 1, f), u.x);
    float x4 = glm::mix(dot_grad(0, 1, 1, f), dot_grad(1, 1, 1, f), u.x);

    float y1 = glm::mix(x1, x2, u.y);
    float y2 = glm::mix(x3, x4, u.y);

    return glm::mix(y1, y2, u.z) * 0.5f + 0.5f;
}

inline float simplex_noise_2d(const glm::vec2& p, uint32_t seed = 0)
{
    const float F2 = 0.366025403784f;
    const float G2 = 0.211324865405f;

    glm::vec2 s = p + (p.x + p.y) * F2;
    glm::vec2 i = glm::floor(s);
    glm::vec2 t = i - (i.x + i.y) * G2;
    glm::vec2 d0 = p - t;
    glm::vec2 d1 =
        d0 - (d0.x > d0.y ? glm::vec2(1.0f, 0.0f) : glm::vec2(0.0f, 1.0f)) + G2;

    glm::vec2 d2 = d0 - glm::vec2(1.0f) + 2.0f * G2;

    uint32_t ix = static_cast<uint32_t>(i.x);
    uint32_t iy = static_cast<uint32_t>(i.y);

    float n0 = 0.0f, n1 = 0.0f, n2 = 0.0f;

    float t0 = 0.5f - glm::dot(d0, d0);
    if (t0 >= 0.0f) {
        t0 *= t0;
        glm::vec2 g = grad2d(hash(ix, iy + seed));
        n0 = t0 * t0 * glm::dot(g, d0);
    }

    float t1 = 0.5f - glm::dot(d1, d1);
    if (t1 >= 0.0f) {
        t1 *= t1;
        glm::vec2 g = d0.x > d0.y ? grad2d(hash(ix + 1, iy + seed))
                                  : grad2d(hash(ix, iy + 1 + seed));
        n1 = t1 * t1 * glm::dot(g, d1);
    }

    float t2 = 0.5f - glm::dot(d2, d2);
    if (t2 >= 0.0f) {
        t2 *= t2;
        glm::vec2 g = grad2d(hash(ix + 1, iy + 1 + seed));
        n2 = t2 * t2 * glm::dot(g, d2);
    }

    return 35.0f * (n0 + n1 + n2) * 0.5f + 0.5f;
}

inline float fbm_2d(
    const glm::vec2& p,
    int octaves = 6,
    float lacunarity = 2.0f,
    float gain = 0.5f,
    uint32_t seed = 0)
{
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float sum = 0.0f;
    float max_value = 0.0f;

    for (int i = 0; i < octaves; ++i) {
        sum += amplitude * perlin_noise_2d(p * frequency, seed + i);
        max_value += amplitude;
        amplitude *= gain;
        frequency *= lacunarity;
    }

    return sum / max_value;
}

inline float fbm_3d(
    const glm::vec3& p,
    int octaves = 6,
    float lacunarity = 2.0f,
    float gain = 0.5f,
    uint32_t seed = 0)
{
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float sum = 0.0f;
    float max_value = 0.0f;

    for (int i = 0; i < octaves; ++i) {
        sum += amplitude * perlin_noise_3d(p * frequency, seed + i);
        max_value += amplitude;
        amplitude *= gain;
        frequency *= lacunarity;
    }

    return sum / max_value;
}

inline float turbulence_2d(
    const glm::vec2& p,
    int octaves = 6,
    float lacunarity = 2.0f,
    float gain = 0.5f,
    uint32_t seed = 0)
{
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float sum = 0.0f;
    float max_value = 0.0f;

    for (int i = 0; i < octaves; ++i) {
        sum += amplitude *
               std::abs(perlin_noise_2d(p * frequency, seed + i) * 2.0f - 1.0f);
        max_value += amplitude;
        amplitude *= gain;
        frequency *= lacunarity;
    }

    return sum / max_value;
}

inline float turbulence_3d(
    const glm::vec3& p,
    int octaves = 6,
    float lacunarity = 2.0f,
    float gain = 0.5f,
    uint32_t seed = 0)
{
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float sum = 0.0f;
    float max_value = 0.0f;

    for (int i = 0; i < octaves; ++i) {
        sum += amplitude *
               std::abs(perlin_noise_3d(p * frequency, seed + i) * 2.0f - 1.0f);
        max_value += amplitude;
        amplitude *= gain;
        frequency *= lacunarity;
    }

    return sum / max_value;
}

inline float voronoi_2d(const glm::vec2& p, uint32_t seed = 0)
{
    glm::vec2 i = glm::floor(p);
    glm::vec2 f = glm::fract(p);

    float min_dist = 1.0f;

    for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
            glm::vec2 neighbor = glm::vec2(x, y);
            glm::vec2 point = neighbor + random_vec2(i + neighbor, seed) - f;
            float dist = glm::length(point);
            min_dist = glm::min(min_dist, dist);
        }
    }

    return glm::clamp(min_dist, 0.0f, 1.0f);
}

inline float voronoi_3d(const glm::vec3& p, uint32_t seed = 0)
{
    glm::vec3 i = glm::floor(p);
    glm::vec3 f = glm::fract(p);

    float min_dist = 1.0f;

    for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
            for (int z = -1; z <= 1; ++z) {
                glm::vec3 neighbor = glm::vec3(x, y, z);
                glm::vec3 point =
                    neighbor + random_vec3(i + neighbor, seed) - f;
                float dist = glm::length(point);
                min_dist = glm::min(min_dist, dist);
            }
        }
    }

    return glm::clamp(min_dist, 0.0f, 1.0f);
}

inline float ridged_multifractal_2d(
    const glm::vec2& p,
    int octaves = 6,
    float lacunarity = 2.0f,
    float gain = 0.5f,
    float offset = 1.0f,
    uint32_t seed = 0)
{
    float frequency = 1.0f;
    float amplitude = 1.0f;
    float sum = 0.0f;
    float weight = 1.0f;

    for (int i = 0; i < octaves; ++i) {
        float signal = perlin_noise_2d(p * frequency, seed + i) * 2.0f - 1.0f;
        signal = std::abs(signal);
        signal = offset - signal;
        signal *= signal;
        signal *= weight;
        weight = glm::clamp(signal * gain, 0.0f, 1.0f);
        sum += signal * amplitude;
        frequency *= lacunarity;
    }

    return sum / 1.5f;
}

inline float ridged_multifractal_3d(
    const glm::vec3& p,
    int octaves = 6,
    float lacunarity = 2.0f,
    float gain = 0.5f,
    float offset = 1.0f,
    uint32_t seed = 0)
{
    float frequency = 1.0f;
    float amplitude = 1.0f;
    float sum = 0.0f;
    float weight = 1.0f;

    for (int i = 0; i < octaves; ++i) {
        float signal = perlin_noise_3d(p * frequency, seed + i) * 2.0f - 1.0f;
        signal = std::abs(signal);
        signal = offset - signal;
        signal *= signal;
        signal *= weight;
        weight = glm::clamp(signal * gain, 0.0f, 1.0f);
        sum += signal * amplitude;
        frequency *= lacunarity;
    }

    return sum / 1.5f;
}

}  // namespace Noise

RUZINO_NAMESPACE_CLOSE_SCOPE
