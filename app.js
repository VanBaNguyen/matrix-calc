"use strict";

/* ----------  tiny complex helpers  --------------------------------- */
const EPS  = 1e-12;                     // arithmetic tolerance
const ZERO = math.complex(0, 0);

const cAbs  = x          => math.abs(x);
const cAdd  = (a, b)     => math.add(a, b);
const cSub  = (a, b)     => math.subtract(a, b);
const cMul  = (a, b)     => math.multiply(a, b);
const cDiv  = (a, b)     => math.divide(a, b);
const cConj = x          => math.conj(x);
const isZero = (x,eps=EPS) => cAbs(x) < eps;

/* nicer string for table display (2 dp, preserves imaginary part) */
function nice(x) {
  if (math.typeOf(x) === "Complex") {
    const re = Math.abs(x.re) < EPS ? 0 : +x.re.toFixed(2);
    const im = Math.abs(x.im) < EPS ? 0 : +x.im.toFixed(2);
    if (im === 0) return `${re}`;
    if (re === 0) return `${im}i`;
    return `${re}${im >= 0 ? "+" : ""}${im}i`;
  }
  return +(+x).toFixed(2);
}

/* ----------  UI helpers  ------------------------------------------- */
function onModeChange() { generateMatrix(); }
function isAxBMode()    { return document.getElementById("mode-axb").checked; }

/* build / rebuild the input grid */
function generateMatrix() {
  const m    = +document.getElementById("rows").value;
  const n    = +document.getElementById("cols").value;
  const form = document.getElementById("matrix-form");
  form.innerHTML = "";

  const totalCols = isAxBMode() ? n + 1 : n;
  form.style.gridTemplateColumns = `repeat(${totalCols}, 70px)`;

  for (let i = 0; i < m; ++i) {
    for (let j = 0; j < totalCols; ++j) {
      const inp        = document.createElement("input");
      inp.type         = "text";          // allows “3+4i”, “-2i”, etc.
      inp.id           = `cell-${i}-${j}`;
      inp.autocomplete = "off";
      inp.value = (j === n)           ? "0" :      // b-vector default
                  (i === j && j < n) ? "1" : "0";  // identity default
      if (isAxBMode() && j === n) inp.classList.add("b-vector-input");
      form.appendChild(inp);
    }
  }
  document.getElementById("calc-btn").style.display = "inline-block";
  document.getElementById("output").innerHTML = "";
}

/* read inputs ⇒ plain JS arrays (numbers or math.Complex) */
function getMatrix() {
  const m   = +document.getElementById("rows").value;
  const n   = +document.getElementById("cols").value;
  const axb = isAxBMode();

  const A = [], b = [];
  for (let i = 0; i < m; ++i) {
    const row = [];
    for (let j = 0; j < (axb ? n + 1 : n); ++j) {
      const str = (document.getElementById(`cell-${i}-${j}`).value || "0").trim();
      let val;
      try { val = math.evaluate(str); } catch { val = 0; }
      if (typeof val === "number") val = +val;        // primitive
      if (axb && j === n) b.push(val);
      else                row.push(val);
    }
    A.push(row);
  }
  return axb ? {A, b} : A;
}

/* pretty-print matrix / vector → HTML */
function toHtml(mat) {
  return `<div class="matrix-display">${
    mat.map(r => "[" + r.map(nice).join(", ") + "]").join("<br>")
  }</div>`;
}

/* ----------  Gauss / REF / RREF  ----------------------------------- */
function ref(mat)  { return _gauss(mat, false); }
function rref(mat) { return _gauss(mat, true ); }

function _gauss(matrix, reduced) {
  const arr  = matrix.map(r => r.slice());           // deep copy
  const rows = arr.length;
  const cols = arr[0].length;
  let lead   = 0;

  for (let r = 0; r < rows && lead < cols; ++r, ++lead) {
    /* 1. find pivot row */
    let i = r;
    while (i < rows && isZero(arr[i][lead])) ++i;
    if (i === rows) { --r; continue; }               // no pivot → next col

    /* 2. swap into place */
    [arr[i], arr[r]] = [arr[r], arr[i]];

    /* 3. scale pivot row to make pivot = 1 */
    const piv = arr[r][lead];
    for (let j = 0; j < cols; ++j) arr[r][j] = cDiv(arr[r][j], piv);

    /* 4. clear column */
    for (let k = 0; k < rows; ++k) {
      if (k === r) continue;
      if (!reduced && k < r) continue;               // REF: only below
      const factor = arr[k][lead];
      if (isZero(factor)) continue;
      for (let j = 0; j < cols; ++j)
        arr[k][j] = cSub(arr[k][j], cMul(factor, arr[r][j]));
    }
  }
  return arr;
}

/* ----------  Solve Ax = b  ---------------------------------------- */
function solveAxEqualsB(A, b) {
  const aug  = A.map((row, i) => row.concat([b[i]]));
  const augR = rref(aug);

  const m = A.length, n = A[0].length;
  const sol = Array(n).fill(ZERO);
  let free  = false;

  /* inconsistency? */
  for (const row of augR) {
    const piv = row.findIndex(x => !isZero(x));
    if (piv === n && !isZero(row[n])) return {type:"none"};
  }

  /* free vars? */
  for (let c = 0; c < n; ++c) {
    const colHas = augR.some(row => !isZero(row[c]));
    if (!colHas) free = true;
  }
  if (free || m < n) return {type:"infinite", rref:augR};

  /* unique solution */
  for (let i = 0; i < n; ++i) sol[i] = augR[i][n];
  return {type:"unique", solution:sol, rref:augR};
}

/* ----------  QR-algorithm (complex)  ------------------------------ */
function dot(u,v)   { return u.reduce((s,ui,idx)=> cAdd(s, cMul(cConj(ui),v[idx])), ZERO); }
function norm(v)    { return Math.sqrt(v.reduce((s,x)=> s + Math.pow(cAbs(x),2), 0)); }

function qrDecomp(A) {
  const n = A.length;
  const R = Array.from({length:n}, () => Array(n).fill(ZERO));
  const Qcols = [];

  const V = [];                              // columns of A
  for (let j = 0; j < n; ++j) V.push(A.map(r => r[j]));

  for (let i = 0; i < n; ++i) {
    let v = V[i].slice();
    for (let j = 0; j < i; ++j) {
      const r = dot(Qcols[j], V[i]);
      R[j][i] = r;
      for (let k = 0; k < n; ++k) v[k] = cSub(v[k], cMul(r, Qcols[j][k]));
    }
    const rii = norm(v);
    R[i][i]   = rii;
    Qcols[i]  = v.map(x => cDiv(x, rii || 1));
  }

  const Q = Array.from({length:n}, (_, r) =>
              Array.from({length:n}, (_, c) => Qcols[c][r]));
  return {Q, R};
}

function matMul(A,B) {
  const n=A.length, m=B[0].length, k=B.length;
  const C = Array.from({length:n}, ()=>Array(m).fill(ZERO));
  for (let i=0;i<n;++i)
    for (let j=0;j<m;++j)
      for (let t=0;t<k;++t)
        C[i][j] = cAdd(C[i][j], cMul(A[i][t], B[t][j]));
  return C;
}

function qrEigenvalues(A, maxIter=150, tol=1e-9) {
  let Ak = A.map(r=>r.slice());
  const n = Ak.length;

  for (let iter=0; iter<maxIter; ++iter) {
    const {Q,R} = qrDecomp(Ak);
    const next  = matMul(R,Q);

    let delta = 0;
    for (let i=0;i<n;++i)
      for (let j=0;j<n;++j)
        delta += Math.abs(cAbs(next[i][j] - Ak[i][j]));
    if (delta < tol) break;
    Ak = next;
  }
  return Ak.map((row,i)=> row[i]);      // diagonal ≈ eigenvalues
}

/* ----------  improved eigenvector finder  ------------------------- */
function findEigenvector(A, λ, tol=1e-6) {
  const n = A.length;

  /* build (A − λI) */
  const M = A.map((row,i) =>
               row.map((v,j)=> cSub(v, (i===j) ? λ : ZERO)));

  /* RREF with relaxed tolerance to expose null-space */
  const R = _gauss(M.map(r=>r.slice()), true);   // own rref

  /* identify pivot columns */
  const pivotCols = new Set();
  for (const row of R) {
    const piv = row.findIndex(x => !isZero(x, tol));
    if (piv !== -1) pivotCols.add(piv);
  }

  /* first free column becomes 1, back-sub for pivots */
  const free = [...Array(n).keys()].filter(j => !pivotCols.has(j));
  if (free.length === 0) return Array(n).fill(ZERO);   // should not happen

  const x = Array(n).fill(ZERO);
  x[free[0]] = math.complex(1,0);

  for (let i = R.length-1; i >= 0; --i) {
    const row = R[i];
    const piv = row.findIndex(v => !isZero(v, tol));
    if (piv === -1) continue;

    let s = ZERO;
    for (let j = piv+1; j < n; ++j)
      s = cAdd(s, cMul(row[j], x[j]));
    x[piv] = cSub(ZERO, s);
  }

  /* normalise for nicer output */
  const vnorm = norm(x);
  return vnorm < EPS ? x : x.map(e => cDiv(e, vnorm));
}

/* ----------  main “Calculate” button ----------------------------- */
function calculateMatrix() {
  const axb = isAxBMode();
  let html  = "";

  if (axb) {
    /* -- solve Ax = b ------------------------------------------- */
    const {A,b} = getMatrix();
    html += `<h3>Matrix A</h3>${toHtml(A)}`;
    html += `<h3>Vector b</h3>${toHtml(b.map(x => [x]))}`;

    const res = solveAxEqualsB(A,b);
    const aug = A.map((r,i)=> r.concat([b[i]]));

    html += `<h3>Augmented REF</h3>${toHtml(ref(aug))}`;
    html += `<h3>Augmented RREF</h3>${toHtml(res.rref || [])}`;

    if (res.type === "unique") {
      html += `<h3>Solution x</h3>${toHtml(res.solution.map(x=>[x]))}`;
    } else if (res.type === "none") {
      html += "<h3>No Solution</h3>";
    } else {
      html += "<h3>Infinite Solutions (free variables)</h3>";
    }

  } else {
    /* -- REF / RREF / eigen ------------------------------------- */
    const M = getMatrix();
    html += `<h3>Input Matrix</h3>${toHtml(M)}`;
    html += `<h3>Row Echelon (REF)</h3>${toHtml(ref(M))}`;
    html += `<h3>Reduced Row Echelon (RREF)</h3>${toHtml(rref(M))}`;

    if (M.length === M[0].length) {      // square ⇒ eigen work
      const eigVals = qrEigenvalues(M);
      html += "<h3>Eigenvalues (approx.)</h3>" +
              eigVals.map((e,i)=>`λ<sub>${i+1}</sub>=${nice(e)}`).join(", ");

      html += "<h3>Eigenvectors (approx.)</h3>";
      eigVals.forEach((λ,i)=>{
        const v = findEigenvector(M, λ);
        html += `λ<sub>${i+1}</sub>: [${v.map(nice).join(", ")}]<br>`;
      });
    }
  }
  document.getElementById("output").innerHTML = html;
}

/* ----------  initial setup  ------------------------------------- */
window.onload = onModeChange;


// === Theme Switcher ===
const THEME_KEY = 'matrixcalc-theme';
const themes = [
  'blue-glass',
  'forest-night',
  'crimson-red',
  'matrix-retro',
  'matrix-green'
];
function applyTheme(themeName) {
  const body = document.body;
  for (let t of themes) body.classList.remove(`theme-${t}`);
  body.classList.add(`theme-${themeName}`);
  // Button visual highlight
  document.querySelectorAll('.theme-btn').forEach(btn => {
    btn.classList.toggle('selected', btn.dataset.theme === themeName);
  });
  localStorage.setItem(THEME_KEY, themeName);
}

// On load, apply saved theme or default
window.addEventListener('DOMContentLoaded', () => {
  const saved = localStorage.getItem(THEME_KEY) || 'blue-glass';
  applyTheme(saved);
  // Setup event listeners for theme buttons
  document.querySelectorAll('.theme-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      applyTheme(btn.dataset.theme);
    });
  });
});
