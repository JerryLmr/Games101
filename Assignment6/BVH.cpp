#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include "BVH.hpp"
#include "Intersection.hpp"

static const char* splitMethodName(BVHAccel::SplitMethod splitMethod)
{
    return splitMethod == BVHAccel::SplitMethod::SAH ? "SAH" : "NAIVE";
}

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    auto start = std::chrono::steady_clock::now();
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    auto stop = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    printf(
        "\rBVH Generation complete (%s): \nTime Taken: %lld ms\n\n",
        splitMethodName(splitMethod), static_cast<long long>(ms));
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        auto getCoord = [](const Vector3f& vector, int dim) {
            if (dim == 0)
                return vector.x;
            if (dim == 1)
                return vector.y;
            return vector.z;
        };

        auto sortByCentroid = [&objects, dim, &getCoord]() {
            std::sort(objects.begin(), objects.end(), [&](auto f1, auto f2) {
                return getCoord(f1->getBounds().Centroid(), dim) <
                       getCoord(f2->getBounds().Centroid(), dim);
            });
        };

        std::vector<Object*> leftshapes;
        std::vector<Object*> rightshapes;

        auto splitByMiddle = [&]() {
            sortByCentroid();

            auto beginning = objects.begin();
            auto middling = objects.begin() + (objects.size() / 2);
            auto ending = objects.end();

            leftshapes = std::vector<Object*>(beginning, middling);
            rightshapes = std::vector<Object*>(middling, ending);
        };

        if (splitMethod == SplitMethod::NAIVE ||
            getCoord(centroidBounds.pMax, dim) ==
                getCoord(centroidBounds.pMin, dim) ||
            bounds.SurfaceArea() <= 0) {
            splitByMiddle();
        }
        else {
            constexpr int nBuckets = 12;
            struct BucketInfo {
                int count = 0;
                Bounds3 bounds;
            };

            std::array<BucketInfo, nBuckets> buckets;

            auto getBucketIndex = [&](Object* object) {
                int bucketIndex =
                    nBuckets *
                    getCoord(centroidBounds.Offset(
                                 object->getBounds().Centroid()),
                             dim);
                if (bucketIndex == nBuckets)
                    bucketIndex = nBuckets - 1;
                if (bucketIndex < 0)
                    bucketIndex = 0;
                if (bucketIndex >= nBuckets)
                    bucketIndex = nBuckets - 1;
                return bucketIndex;
            };

            for (auto object : objects) {
                int bucketIndex = getBucketIndex(object);
                assert(bucketIndex >= 0 && bucketIndex < nBuckets);
                buckets[bucketIndex].count++;
                buckets[bucketIndex].bounds =
                    Union(buckets[bucketIndex].bounds, object->getBounds());
            }

            std::array<float, nBuckets - 1> cost;
            for (int i = 0; i < nBuckets - 1; ++i) {
                Bounds3 leftBounds, rightBounds;
                int leftCount = 0, rightCount = 0;

                for (int j = 0; j <= i; ++j) {
                    if (buckets[j].count == 0)
                        continue;
                    leftBounds = Union(leftBounds, buckets[j].bounds);
                    leftCount += buckets[j].count;
                }
                for (int j = i + 1; j < nBuckets; ++j) {
                    if (buckets[j].count == 0)
                        continue;
                    rightBounds = Union(rightBounds, buckets[j].bounds);
                    rightCount += buckets[j].count;
                }

                if (leftCount == 0 || rightCount == 0) {
                    cost[i] = std::numeric_limits<float>::infinity();
                }
                else {
                    cost[i] =
                        1.0f +
                        (leftCount * leftBounds.SurfaceArea() +
                         rightCount * rightBounds.SurfaceArea()) /
                            bounds.SurfaceArea();
                }
            }

            float minCost = cost[0];
            int minCostSplitBucket = 0;
            for (int i = 1; i < nBuckets - 1; ++i) {
                if (cost[i] < minCost) {
                    minCost = cost[i];
                    minCostSplitBucket = i;
                }
            }

            if (!std::isfinite(minCost)) {
                splitByMiddle();
            }
            else {
                auto middling =
                    std::partition(objects.begin(), objects.end(),
                                   [&](Object* object) {
                                       return getBucketIndex(object) <=
                                              minCostSplitBucket;
                                   });

                if (middling == objects.begin() || middling == objects.end()) {
                    splitByMiddle();
                }
                else {
                    leftshapes = std::vector<Object*>(objects.begin(), middling);
                    rightshapes = std::vector<Object*>(middling, objects.end());
                }
            }
        }

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    Intersection inter;
    std::array<int,3> dirIsNeg = {
        ray.direction.x < 0,
        ray.direction.y < 0,
        ray.direction.z < 0
    };
    if (!node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg))
        return inter;
    if (node->left == nullptr && node->right == nullptr) 
        return node->object->getIntersection(ray);

    Intersection Left = getIntersection(node->left, ray);
    Intersection Right = getIntersection(node->right, ray);

    if (Left.happened && Right.happened) 
        return Left.distance < Right.distance ? Left : Right;

    if (Left.happened) 
        return Left;

    return Right;
}
