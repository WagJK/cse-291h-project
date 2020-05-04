#version 330 core

uniform vec3 fragColor;
out vec4 FragColor;

void main()
{
	// linearly interpolate between both textures (80% container, 20% awesomeface)
	FragColor = vec4(fragColor, 1.f);
}