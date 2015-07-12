#version 120   
#extension EXT_gpu_shader4 : require    

varying vec3 normal, lightDir[3], eyeVec;

void main() 
{  
	float fScale = 0.5;
    vec4 final_color = vec4(0, 0, 0, 0);
    
    for (int i = 0; i < 3; i ++)
    {
		final_color += gl_FrontLightModelProduct.sceneColor * gl_FrontMaterial.ambient * fScale + 
					   gl_LightSource[i].ambient * gl_FrontMaterial.ambient * fScale;
							
		vec3 N = normalize(normal);
		vec3 L = normalize(lightDir[i]);
	
		float lambertTerm = dot(N,L);
	
		if(lambertTerm > 0.0)
		{
			final_color += gl_LightSource[i].diffuse * 
				           gl_FrontMaterial.diffuse * 
						   lambertTerm * fScale;	
		
			vec3 E = vec3(0, 0, 1);
			vec3 R = reflect(-L, N);
			float specular = pow( max(dot(R, E), 0.0), 30 );
			final_color += gl_LightSource[i].specular * 
				           vec4(1, 1, 1, 1) * 
						   specular;	
		}
    }
	
	gl_FragColor = final_color;	
}