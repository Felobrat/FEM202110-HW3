The new implementation of the nodal stress averaging we implemeted a four side square like fillipa chp. 28 and we calculated the tributary area for each quad4 node ande the four element each share the node
inclugind the results of the 3 types of mesh, a fine, medium and coarse. 
![deformada-gruesa](https://user-images.githubusercontent.com/53713268/118558973-4f81fa00-b735-11eb-97d9-a20d0f4b0bf3.png)
"The coarse mesh from python"
![deformada-media](https://user-images.githubusercontent.com/53713268/118559099-79d3b780-b735-11eb-8581-913ac16b9265.png)
"The medium mesh from python"
![deformada-fina](https://user-images.githubusercontent.com/53713268/118559218-9b34a380-b735-11eb-971b-c372ea78ce3c.png)
"The fine from python"\
The displacement of each mesh was:

- Coarse:
![desplazamiento_coarse](https://user-images.githubusercontent.com/53713268/118560476-8bb65a00-b737-11eb-979b-aa247c4f9eba.JPG)

- Medium:
![desplazamiento](https://user-images.githubusercontent.com/53713268/118559717-56f5d300-b736-11eb-849c-87ae766fb1ee.JPG)

-Fine:
[desplazamiento_fino](https://user-images.githubusercontent.com/53713268/118560661-d89a3080-b737-11eb-8c87-4a8c48dfcd8d.JPG)

the stresses for sigmax was:
- Coarse:
![Sigma_x_grueso](https://user-images.githubusercontent.com/53713268/118569085-28ccbf00-b747-11eb-9c1f-9ab417577826.jpeg)


- Medium:
![sigma_x_medio](https://user-images.githubusercontent.com/53713268/118569138-4bf76e80-b747-11eb-9afe-84a5c9286e38.jpeg)


- Fine:
![image](https://user-images.githubusercontent.com/53713268/118561323-d2f11a80-b738-11eb-8aed-439d54289325.png)


For the nodal stress averingig we take a subsquare to calculated the stress in each point of gauss (4 for element) with this 4 point, we can calculated the influence of (triutary area) of each element wich use te node and influence it.
 
