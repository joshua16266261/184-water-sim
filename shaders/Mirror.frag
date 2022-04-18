#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
	vec3 v_pos_3 = vec3(v_position.x, v_position.y, v_position.z);
	vec3 v_n_3 = vec3(v_normal.x, v_normal.y, v_normal.z);
	vec3 wo = v_pos_3 - u_cam_pos;
	vec3 wi = wo - 2.0 * dot(wo, v_n_3) * normalize(v_n_3);
	out_color = texture(u_texture_cubemap, wi);
}
