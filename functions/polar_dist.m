function delta_sigma = polar_dist(theta_1, theta_2, phi_1, phi_2)
   delta_theta = abs(theta_1 - theta_2); 
   delta_phi= abs(phi_1 - phi_2); 

   cos_phi_1 = cos(phi_1);
   cos_phi_2 = cos(phi_2);
   sin_phi_1 = sin(phi_1);
   sin_phi_2 = sin(phi_2);

   cos_delta_theta = cos(delta_theta);
   
   a = (cos_phi_2 * sin(delta_theta))^2;
   b = (cos_phi_1 * sin_phi_2 - sin_phi_1 * cos_phi_2 * cos_delta_theta)^2;
   c = sin_phi_1 * sin_phi_2 + cos_phi_1 * cos_phi_2 * cos_delta_theta;

   delta_sigma = atan2( sqrt(a + b), c);
end
