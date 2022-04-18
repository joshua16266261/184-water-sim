#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
	return texture(u_texture_2, uv).r;
}

void main() {
  // YOUR CODE HERE
  
	vec3 v_normal_3 = vec3(v_normal.x, v_normal.y, v_normal.z);
	vec3 v_tangent_3 = vec3(v_tangent.x, v_tangent.y, v_tangent.z);
	vec3 b = cross(v_normal_3, v_tangent_3);
	mat3 tbn = mat3(v_tangent_3, b, v_normal_3);
	
	float du = (h(vec2(v_uv.x + 1.0 / u_texture_2_size.x, v_uv.y)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
	float dv = (h(vec2(v_uv.x, v_uv.y + 1.0 / u_texture_2_size.y)) - h(v_uv)) * u_height_scaling * u_normal_scaling;
	
	vec3 n0 = vec3(-du, -dv, 1.0);
	vec3 nd = tbn * n0;
	
	float p = 69.0;
	float ka = 0.1;
	float kd = 0.9;
	float ks = 0.9;
	
	vec3 Ia = vec3(1.0, 1.0, 1.0);
	
	vec3 out_3 = ka * Ia;
	
	vec3 v_pos_3 = vec3(v_position.x, v_position.y, v_position.z);
	vec3 v_n_3 = normalize(nd);
	
	vec3 v = normalize(u_cam_pos - v_pos_3);
	vec3 l = normalize(u_light_pos - v_pos_3);
	vec3 h = normalize(v + l);
	
	
  
	// Diffuse
	out_3 += kd * u_light_intensity / pow(distance(v_pos_3, u_light_pos), 2.0) * max(0.0, dot(v_n_3, l));
	
	
	// Specular
	out_3 += ks * u_light_intensity / pow(distance(v_pos_3, u_light_pos), 2.0) * pow(max(0.0, dot(v_n_3, h)), p);
	
	
	out_color = vec4(out_3.x, out_3.y, out_3.z, 1.0);
}

