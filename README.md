# Computer-Graphics-Engine
This repository represents a Computer Graphics Engine done from scratch entirely in C++. 

# Showcase 
Here you can see a couple of pictures generated using only our Engine:
![image](https://user-images.githubusercontent.com/99887781/220395055-7f14e49e-6923-4ddc-b4f6-0cd398d1b4de.png)
![image](https://user-images.githubusercontent.com/99887781/220393876-e813252d-f736-4db1-bad6-bd50212ff5d3.png)
![image](https://user-images.githubusercontent.com/99887781/220395087-996cc3a4-d1f8-495f-bff8-ef6db24a8c4c.png)

![image](https://user-images.githubusercontent.com/99887781/220394320-35202ed4-ea51-48cf-a2dc-6339a8448df7.png)
![image](https://user-images.githubusercontent.com/99887781/220393944-8479673c-774a-43eb-975e-ecf5de97a6eb.png)

# Features
The features implemented are as follows:

## Shading using Phong shading model
### Feature:
![image](https://user-images.githubusercontent.com/99887781/220385121-f0352e52-bbaa-465a-9a73-2a5b5e6645da.png)
### Debugger:
Here a ray is shot at the right wall of the box and it is colored green, as the color of the intersection point
![image](https://user-images.githubusercontent.com/99887781/220385433-9baf8a60-721b-4c8b-b506-7a1d0eaf2a67.png)
![image](https://user-images.githubusercontent.com/99887781/220385689-4442d6b8-f9cb-40b9-84e8-54aa85687516.png)

## Area lights

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220385907-462c5174-01e5-41ec-8f5f-dc236c8035dc.png)
### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220386014-54bb0546-cfdc-4354-aa49-6ccf1dc483b4.png)
![image](https://user-images.githubusercontent.com/99887781/220386093-0d86bfb7-b52a-4771-9623-fc97cc0e8992.png)

## Normal interpolation with barycentric coordinates

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220387655-18a9975c-9c10-4ea9-99c7-1a2fa59ab636.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220387692-f9581f01-5e53-49c5-b4ec-c72a74e82006.png)

## Texture mapping

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220387781-3ea844ff-cdbe-4053-95c6-abf423d8cb64.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220387857-f9316092-24aa-49ce-acf8-df3080e93ea6.png)

## Recursive Ray-Tracer

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220387963-583ebd22-ede0-4a34-a940-f2fe16d41b31.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220388019-57fd41a9-3e96-46cb-99ea-0996401e21c9.png)

## Hard shadows

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220388084-ceea62cc-d10d-4041-975c-933e8df9eddb.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220388153-55682af3-36ce-4763-88e6-38367b83ddaf.png)

## Bounding Volume Hierarchy - Creation

### Implementation:

The BVH is implemented fully inside the bounding_volume_hierarchy.cpp and the
corresponding header file.
The creation goes the following way:
1. For each triangle and sphere in the scene create a Primitive struct. This structure
stores the primitive’s centroid as well as the data of the primitive itself. For the
triangle it collects the triangle’s mesh index and the indices of the vertices. For the
sphere it saves the sphere’s index in the scene.
2. Now we call the createBVH function. It is recursively defined to operate on a range of
primitives inside the primitives vector.
3. If the range contains only one primitive or we have achieved the maximum depth
(defined to be 16 levels), we are creating the leaf node. The node struct stores an
axis aligned bounding box that holds all the primitives inside. Additionally it has some
more information - a boolean marker for leaf nodes, the range into the primitives
vector (of the contained primitives) and for non-leaf nodes the indices of the left and
right child nodes. Once we have all the information about the node, we save it in the
nodes vector and return the index we stored it at.
4. For the recursive case we call the function splitFunc that for the standard feature
works the following way. We pass the primitives vector, the range into this vector we
want to operate on and the depth we’re at. For the standard BVH the splitFunc is set
to be splitStandard function. Based on the depth value we decide the axis to perform
the split on. Then we find the median primitive in the range by doing a partial sort
(std::nth_element). Now we have all the “smaller” primitives to the left of the median
and all the “greater” to the right. Lastly we return the index on which we performed
the split.
5. Back in createBVH we call the function recursively on the left and right ranges. We
get back the indices where the child nodes got stored at. Now we have all the
required information to save our internal node, so we do it. Lastly we return the index
the node was stored at.
6. Now the bounding volume hierarchy is finished. We save the root node’s index to be
used for traversal.

### Performance:
![image](https://user-images.githubusercontent.com/99887781/220388473-972857a6-2818-402a-90f2-8b62acea4c66.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220388574-db430593-0975-4710-b93f-e0c801c88744.png)

## Bounding Volume Hierarchy - Traversal

### Implementation:
In the method intersect() in bounding_volume_hierarchy.cpp a stack is being used for the
traversal of the hierarchy, as it will contain all the intersected aabb’s found in the ray’s path.
We start with the root node that is first pushed as the first element in the stack. We have a
while loop that iterates until the stack becomes empty(no more unvisited aabb’s). Every
iteration we pop the top element from the stack, check whether that popped node is a leaf
node, if it is a leaf node we iterate over all the primitives that leaf node contains and use the
intersect method provided by the library to check whether an intersection has happened. If a
node is not a leaf node, then we take its two children (we use indexes here into the
std::vector<Node> nodes list to retrieve the corresponding children). After getting the
children we check if the children are being intersected by the upcoming ray. If they are, we
push them to the stack. Traversal checks all the nodes and updates the closest primitive so
far hit. We have defined a Primitive to be a Triangle or a Sphere, and thus the traversal
handles both cases very nicely. Optimization has been achieved through the use of indexes
and const references, which has improved the speed of rendering images. In the
getIntersecting() method we pass the iterators for the first primitive in the box and the last
primitive in the box, and by using modern c++ we iterate over all of them and update the
optional primitive in case of a hit.

### Performance:
![image](https://user-images.githubusercontent.com/99887781/220389036-98c236eb-29d7-42a1-8b43-07b780d28eca.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220389161-55d52920-efbc-4373-b06d-a6f4ec367775.png)

## SAH + Binning

### Feature and Debugger:
![image](https://user-images.githubusercontent.com/99887781/220389356-12a5d40a-5813-4f3a-b346-30694b08f7db.png)

## Bloom filter on the final image

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220389876-41aa88d6-2af5-481e-a4d6-7c888600a454.png)

## Casting multiple rays per pixel

### Feature:
<img width="497" alt="image" src="https://user-images.githubusercontent.com/99887781/220390067-aab7da54-30e2-4360-9ab4-d6966ec758f0.png">
![image](https://user-images.githubusercontent.com/99887781/220390109-a080919e-7df9-40ca-a817-9d311ba2bdcf.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220390213-57e8c481-3f55-4ff2-a8d7-52a772850f78.png)

## Bilinear interpolation

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220390276-49401460-f611-4685-bdf0-0c376c22d96b.png)
![image](https://user-images.githubusercontent.com/99887781/220390317-74d5f75d-753d-4e8c-8dae-cccd18a5f5da.png)

## Transparency

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220390403-25bdfa63-9257-4f92-8e36-86adcd1dcad1.png)
![image](https://user-images.githubusercontent.com/99887781/220390463-f6bac639-4262-4565-8489-1e3458b975b2.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220390497-e15de888-ed4a-47e3-870d-dc11d3bbde17.png)

## Glossy reflections

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220390568-dffba0eb-551c-473a-bf31-5e6b11676de6.png)
![image](https://user-images.githubusercontent.com/99887781/220390612-23095377-41be-4792-81ad-4561f2e77903.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220390709-b11552e2-2765-496b-9758-efc877d1539d.png)

## Depth of field

### Feature:
![image](https://user-images.githubusercontent.com/99887781/220390760-7cb41417-fc31-4d6e-bd91-c8c6959da4d2.png)
![image](https://user-images.githubusercontent.com/99887781/220390783-842f2386-fbe4-45ab-86f8-6e051c79943e.png)

### Debugger:
![image](https://user-images.githubusercontent.com/99887781/220390875-83017e4e-7eef-4e9e-a3e4-2c71cd49b033.png)
![image](https://user-images.githubusercontent.com/99887781/220390894-2768dd02-24ca-4958-b217-0ba802974dc2.png)




