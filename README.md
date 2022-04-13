# ComplexImageMethod
<h3>Image method for shoebox rooms with spherical wave scattering</h3>

<p>It has been shown by <a href="https://www.jstage.jst.go.jp/article/ast/26/2/26_2_145/_pdf">Lam</a> that the image method with spherical wave reflection coefficients outperforms the image method with 
plane wave reflection coefficients when the surface impedance is far from rigid and the angle of incidence is far from normal. It has been shown to have the same 
accuracy as the boundary element method, when the room is rectangular with uniform admittance. However, this does not hold true in rooms with
complicated geometries nor model scattering due to surface roughness.</p>

<p> This implementation takes in as input room dimensions, source position, microphone array positions, wall impedance and wave numbers to evaluate the solution on,
  and returns as output the total pressure field at the mic array locations.<p>
