#include "Marmot/DisplacementFBarParticleSQCNIxSDI.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool DisplacementFBarParticleSQCNIxSDI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarSQCNIxSDI/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarParticleSQCNIxSDI<
                          2,
                          4 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementFBarParticleSQCNIxSDI< 2, 4 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementFBarParticleSQCNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarSQCNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarParticleSQCNIxSDI<
                          3,
                          8 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementFBarParticleSQCNIxSDI< 3, 8 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementFBarParticleRS_SNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarRS-SNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementFBarParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::
                                                                              RotationAndPrincipalStretch );
                      } );

  const static bool DisplacementFBarParticleR_SNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarR-SNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementFBarParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::RotationOnly );
                      } );

  const static bool DisplacementFBarParticleSNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarSNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementFBarParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::None );
                      } );

} // namespace Marmot::Meshfree
