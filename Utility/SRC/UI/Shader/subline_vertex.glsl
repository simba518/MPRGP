#version 120   
#extension EXT_gpu_shader4 : require    


varying vec3 normal, lightDir[3];

void main() 
{  
    normal = gl_NormalMatrix * gl_Normal;

	vec3 vVertex = vec3(gl_ModelViewMatrix * gl_Vertex);

	for (int i = 0; i < 3; i ++)
	{
		lightDir[i] = vec3(gl_LightSource[i].position.xyz - vVertex);
	}

	gl_Position = gl_ModelViewProjectionMatrix * (gl_Vertex + vec4(gl_Normal * 0.004, 0.0));  	
}