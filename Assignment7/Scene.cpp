//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include "Intersection.hpp"
#include "Vector.hpp"
#include "global.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TODO Implement Path Tracing Algorithm here
    Intersection inter = intersect(ray); // p = inter.coords
    Vector3f N = inter.normal;
    if (!inter.happened)
        return Vector3f(0.0f);

    if (inter.m->hasEmission())
        return inter.m->getEmission();
    
    Vector3f wo = ray.direction;

    Intersection lightInter;
    float pdf_light = 0.0f;
    sampleLight(lightInter, pdf_light);

    Vector3f x = lightInter.coords; //光源随机点x
    Vector3f toLight = x - inter.coords;
    float lightDistance = toLight.norm();
    Vector3f ws = normalize(toLight);
    Vector3f NN = lightInter.normal;
    Vector3f emit = lightInter.emit;

    Ray pTox(inter.coords, ws);
    Intersection shadowInter = intersect(pTox);

    Vector3f L_dir = 0.0f;
    if (shadowInter.happened && fabs(lightDistance - shadowInter.distance)< EPSILON) {
        float cosSurface = std::max(0.0f, dotProduct(ws, N));
        float cosLight = std::max(0.0f, dotProduct(-ws, NN));
        L_dir = emit * inter.m->eval(wo, ws, N) * cosSurface * cosLight 
                      / (lightDistance * lightDistance) / pdf_light;
    }

    Vector3f L_indir = 0.0f;
    if (get_random_float() < RussianRoulette) {
        Vector3f wi = normalize(inter.m->sample(wo, N));
        float pdf = inter.m->pdf(wo, wi, N);

        Ray bounceRay(inter.coords, wi);
        Intersection nextInter = intersect(bounceRay);
        if (nextInter.happened && !nextInter.m->hasEmission() && pdf > 0.0f) {
            float cosine = std::max(0.0f, dotProduct(wi, N));

            L_indir = castRay(bounceRay, depth + 1) * inter.m->eval(wo,wi,N) 
                     * cosine / pdf / RussianRoulette;

        }
    }
    return L_dir + L_indir;
}