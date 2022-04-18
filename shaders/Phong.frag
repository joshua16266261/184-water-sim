#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
	float p = 69.0;
	float ka = 0.1;
	float kd = 0.9;
	float ks = 0.9;
	
	vec3 Ia = vec3(1.0, 1.0, 1.0);
	
	vec3 out_3 = ka * Ia;
	
	vec3 v_pos_3 = vec3(v_position.x, v_position.y, v_position.z);
	vec3 v_n_3 = normalize(vec3(v_normal.x, v_normal.y, v_normal.z));
	
	vec3 v = normalize(u_cam_pos - v_pos_3);
	vec3 l = normalize(u_light_pos - v_pos_3);
	vec3 h = normalize(v + l);
	
	
  
	// Diffuse
	out_3 += kd * u_light_intensity / pow(distance(v_pos_3, u_light_pos), 2.0) * max(0.0, dot(v_n_3, l));
	
	
	// Specular
	out_3 += ks * u_light_intensity / pow(distance(v_pos_3, u_light_pos), 2.0) * pow(max(0.0, dot(v_n_3, h)), p);
	
	
	out_color = vec4(out_3.x, out_3.y, out_3.z, 1.0);
}

