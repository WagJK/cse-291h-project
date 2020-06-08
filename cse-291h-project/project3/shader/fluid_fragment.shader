#version 330 core

uniform vec3 fragColor;
out vec4 FragColor;

const float PI = 3.14159265;

void main()
{
	// linearly interpolate between both textures (80% container, 20% awesomeface)
	
	vec3 lightDir = vec3(0.3, 0.3, 0.9);
	vec3 N;
	N.xy = gl_PointCoord * 2.0 - vec2(1.0);

	float mag = dot(N.xy, N.xy);
	if (mag > 1.0) discard;
	N.z = sqrt(1.0 - mag);

	float diffuse = max(0.0, dot(lightDir, N));
	FragColor = vec4(fragColor, 1.f) * diffuse;
}