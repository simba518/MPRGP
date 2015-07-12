//--------------------------------------------------------------------------------------
// Order Independent Transparency with Dual Depth Peeling
//
// Author: Louis Bavoil
// Email: sdkfeedback@nvidia.com
//
// Copyright (c) NVIDIA Corporation. All rights reserved.
//--------------------------------------------------------------------------------------

uniform samplerRECT DepthTex;

void main(void)
{
	// Bit-exact comparison between FP32 z-buffer and fragment depth
	float frontDepth = textureRect(DepthTex, gl_FragCoord.xy).r;
	if (gl_FragCoord.z <= frontDepth) {
		discard;
	}
	
	gl_FragColor = gl_Color;
}
