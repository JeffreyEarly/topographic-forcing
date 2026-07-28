function audit = auditBoundaryDynamicalDescriptor(self)
% Audit the pressure-retaining boundary-dynamical descriptor.
%
% The construction uses primitive hydrostatic F--G coordinates internally,
% retains pressure through the complete saddle-point pencil, and maps the
% finite prognostic modes back to the unchanged public Galerkin layout.
% Milestone 5 supports flat and spatially uniform terrain.
%
% ```matlab
% audit = problem.auditBoundaryDynamicalDescriptor();
% ```
%
% - Topic: Audit boundary-dynamical evolution
% - Declaration: audit = auditBoundaryDynamicalDescriptor()
% - Returns audit: descriptor matrices, reduced modes, and compatibility diagnostics
audit = buildBoundaryDynamicalDescriptorAudit(self);
end
