#version 330 core
layout(location = 0) in vec3 aPos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float pointScale;
uniform float pointRadius;

void main()
{	vec3 posEye = vec3(view * model * vec4(aPos, 1.0));
	float dist = length(posEye);

	gl_PointSize = pointRadius * (pointScale / dist);
	gl_Position = projection * view * model * vec4(aPos, 1.0f);
}