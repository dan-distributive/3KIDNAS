'use strict';

// =============================================================================
// GeometryEstimates.js
// Faithful port of FitDriverScripts/GeometryEstimates.py -- computes the
// Inc/PA estimates 3KIDNAS derives from a SoFiA catalogue entry. Pure
// arithmetic, no Fortran/native dependency, so this is a direct line-by-line
// port, not a re-derivation.
//
// This project's actual usage (SoFiA_Driver.py:LoadSoFiAOutput) hardcodes
// IncMethod='WALLABY_Like' and passes the SAME catalogue PA value in for
// both kinPA and ellPA (both are sourced from the catalogue's kin_pa
// column -- see ParseSoFiACatalog.js's own header comment for why), so
// PAEstimate's nan-fallback branch is effectively dead in practice here.
// Ported faithfully anyway since it's cheap and this file may see other
// IncMethod/PA inputs later.
// =============================================================================

/**
 * @param {{ell_maj: number, ell_min: number, kin_pa: number, ell_pa: number}} sofiaCat
 * @param {number} beamPix - beam size in pixels
 * @param {'UnCorrected'|'BeamCorrected'|'WALLABY_Like'} incMethod
 * @returns {{paEst: number, incEst: number}}
 */
function getGeometryEstimates(sofiaCat, beamPix, incMethod) {
  const paEst = paEstimate(sofiaCat.kin_pa, sofiaCat.ell_pa);
  const incEst = incEstimate(sofiaCat, beamPix, incMethod);
  return { paEst, incEst };
}

function beamCorrSize2(size, beam) {
  return size * size - beam * beam;
}

function incEstimate(sofiaCat, beamPix, incMethod) {
  let minUse2, majUse2;
  if (incMethod === 'UnCorrected') {
    minUse2 = sofiaCat.ell_min ** 2;
    majUse2 = sofiaCat.ell_maj ** 2;
  } else if (incMethod === 'BeamCorrected') {
    minUse2 = beamCorrSize2(sofiaCat.ell_min, beamPix);
    majUse2 = beamCorrSize2(sofiaCat.ell_maj, beamPix);
  } else if (incMethod === 'WALLABY_Like') {
    minUse2 = beamCorrSize2(2 * sofiaCat.ell_min, beamPix);
    majUse2 = beamCorrSize2(2 * sofiaCat.ell_maj, beamPix);
  } else {
    throw new Error('No valid inclination estimate method provided.');
  }

  if (majUse2 < 0) majUse2 = beamPix ** 2;

  let incEst;
  if (minUse2 <= 0) {
    incEst = 89;
  } else {
    incEst = (Math.acos(Math.sqrt(minUse2) / Math.sqrt(majUse2)) * 180) / Math.PI;
  }
  return incEst;
}

function paEstimate(kinPA, ellPA) {
  return Number.isNaN(kinPA) ? ellPA : kinPA;
}

// BOTH assignments, unconditionally -- see ParseSoFiACatalog.js's identical
// comment for the full reasoning (two real, independent consumers: every
// require()-based context needs module.exports set unconditionally;
// galaxy-fit.html's raw-<script>-tag local-estimate loader needs
// globalThis.GeometryEstimates, which the previous if/else version of this
// code dropped whichever branch didn't fire).
if (typeof module !== 'undefined') {
  module.exports = { getGeometryEstimates, incEstimate, paEstimate, beamCorrSize2 };
}
if (typeof globalThis !== 'undefined') {
  globalThis.GeometryEstimates = { getGeometryEstimates, incEstimate, paEstimate, beamCorrSize2 };
}
