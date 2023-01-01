//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


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

// // Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Vector3f dirLight(0.0f, 0.0f, 0.0f), indirLight(0.0f, 0.0f, 0.0f);
    
    // get intersection p
    Intersection pInter = intersect(ray);

    // if hit an object
    if(pInter.happened)
    {
        // hit the light
        if(pInter.m->hasEmission())
        {
            return pInter.m->getEmission();
        }

        Vector3f p = pInter.coords;
        Vector3f wo = -ray.direction;
        Vector3f N = pInter.normal;
        
        // direct light
        Intersection lightInter;
        float pdf;
        sampleLight(lightInter, pdf);
        Vector3f x = lightInter.coords;
        Vector3f NN = lightInter.normal;
        Vector3f ws = normalize(x-p);
        Vector3f emit = lightInter.emit;
        
        // check if the light is blocked
        Ray blockRay(p, ws);
        Intersection blockInter = intersect(blockRay);
        
        if((blockInter.coords - x).norm() < 0.01f)
        {
            auto f_r = pInter.m->eval(wo, ws, N);
            auto cos = dotProduct(ws, N);
            auto cos_x = dotProduct(-ws, NN);
            auto xp2 = dotProduct(x - p, x - p);
            dirLight = emit * f_r * cos * cos_x / xp2 / pdf;
        }
        
        // indirect light
        if(get_random_float() < RussianRoulette)
        {
            Vector3f wi = normalize(pInter.m->sample(wo, N));
            
            Ray indirRay(p, wi);
            Intersection indirInter = intersect(indirRay);
            
            if(indirInter.happened && !indirInter.m->hasEmission())
            {
                indirLight = castRay(Ray(p, wi), depth) * pInter.m->eval(wo, wi, N) * dotProduct(wi, N) / pInter.m->pdf(wo, wi, N) / RussianRoulette;
            }
        }
    }
    
    return dirLight + indirLight;
}