#version 330

// The camera's position in world-space
uniform vec3 u_cam_pos;

// Color
uniform vec4 u_color;

// Properties of the single point light
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

// We also get the uniform texture we want to use.
uniform sampler2D u_texture_1;

// These are the inputs which are the outputs of the vertex shader.
in vec4 v_position;
in vec4 v_normal;

// This is where the final pixel color is output.
// Here, we are only interested in the first 3 dimensions (xyz).
// The 4th entry in this vector is for "alpha blending" which we
// do not require you to know about. For now, just set the alpha
// to 1.
out vec4 out_color;

void main() {
  // YOUR CODE HERE

	vec3 v_pos_3 = vec3(v_position.x, v_position.y, v_position.z);
	vec3 v_n_3 = vec3(v_normal.x, v_normal.y, v_normal.z);
	vec3 out_3 = u_light_intensity / pow(distance(v_pos_3, u_light_pos), 2.0) * max(0.0, dot(normalize(v_n_3), normalize(u_light_pos - v_pos_3)));
	out_color = vec4(out_3.x, out_3.y, out_3.z, 1.0);
}
