function onModeChange() {
  generateMatrix();
}

function isAxBMode() {
  return document.getElementById('mode-axb').checked;
}

function generateMatrix() {
  const m = parseInt(document.getElementById('rows').value);
  const n = parseInt(document.getElementById('cols').value);
  const form = document.getElementById('matrix-form');
  form.innerHTML = '';
  let numCols = isAxBMode() ? n + 1 : n; // add extra column for b

  form.style.gridTemplateColumns = `repeat(${numCols}, 70px)`;

  for (let i = 0; i < m; ++i) {
    for (let j = 0; j < numCols; ++j) {
      const inp = document.createElement('input');
      inp.type = 'number';
      inp.step = 'any';
      inp.id = `cell-${i}-${j}`;
      inp.inputMode = "numeric";
      inp.pattern = "\\d*";
      // Fill diagonal with 1, b (last col) with 0
      inp.value = (j === n) ? 0 : (i === j ? 1 : 0);

      // Add special class if it's the b column in Ax=b mode
      if (isAxBMode() && j === n) {
        inp.classList.add('b-vector-input');
      }

      form.appendChild(inp);
    }
  }
  document.getElementById('calc-btn').style.display = 'inline-block';
  document.getElementById('output').innerHTML = '';
}

function getMatrix() {
  const m = parseInt(document.getElementById('rows').value);
  const n = parseInt(document.getElementById('cols').value);
  let matrix = [];
  let b = [];
  let isAxB = isAxBMode();

  for (let i = 0; i < m; ++i) {
    let row = [];
    for (let j = 0; j < (isAxB ? n + 1 : n); ++j) {
      let val = parseFloat(document.getElementById(`cell-${i}-${j}`).value);
      if (isAxB && j === n) b.push(isNaN(val) ? 0 : val);
      else row.push(isNaN(val) ? 0 : val);
    }
    if (isAxB) matrix.push(row); // matrix is A, b is separate
    else matrix.push(row);       // matrix is full
  }
  if (isAxB) return {A: matrix, b: b};
  else return matrix;
}

function solveAxEqualsB(A, b) {
  // Form the augmented matrix [A | b]
  let aug = A.map((row, i) => row.concat([b[i]]));
  let rrefAug = rref(aug);

  // Check for no solution or infinite solutions
  let m = A.length, n = A[0].length;
  let solution = Array(n).fill(0);
  let freeVars = [];

  for (let i = 0; i < m; ++i) {
    let pivotCol = rrefAug[i].findIndex(val => Math.abs(val) > 1e-8);
    if (pivotCol === n) {
      // Row: [0, 0, ..., 0, c] (c ≠ 0), so no solution
      if (Math.abs(rrefAug[i][n]) > 1e-8) return {type: "none"};
    }
  }

  // If n > m, system is underdetermined, might be infinite
  // For simplicity, just solve unique solution or give "infinite"
  let hasFreeVars = false;
  for (let i = 0; i < n; ++i) {
    let col = rrefAug.map(row => row[i]);
    let nonzeroRows = col.reduce((c, v) => c + (Math.abs(v) > 1e-8 ? 1 : 0), 0);
    if (nonzeroRows === 0) hasFreeVars = true;
  }

  if (hasFreeVars || m < n) return {type: "infinite", rref: rrefAug};

  // Otherwise, extract solution from RREF
  for (let i = 0; i < n; ++i) {
    // Each row corresponds to x_i
    solution[i] = rrefAug[i][n];
  }
  return {type: "unique", solution: solution, rref: rrefAug};
}


function matrixToHtml(matrix) {
  return `<div class="matrix-display">${matrix.map(row =>
    '[' + row.map(x => Number(x).toFixed(2)).join(', ') + ']'
  ).join('<br>')}</div>`;
}

function ref(matrix) {
  matrix = matrix.map(row => row.slice());
  let lead = 0;
  let rowCount = matrix.length;
  let colCount = matrix[0].length;
  for (let r = 0; r < rowCount; r++) {
    if (colCount <= lead) return matrix;
    let i = r;
    while (matrix[i][lead] === 0) {
      i++;
      if (i === rowCount) {
        i = r;
        lead++;
        if (colCount === lead) return matrix;
      }
    }
    let tmp = matrix[i];
    matrix[i] = matrix[r];
    matrix[r] = tmp;

    let val = matrix[r][lead];
    if (val !== 0) {
      for (let j = 0; j < colCount; j++) matrix[r][j] /= val;
    }

    for (let i = r + 1; i < rowCount; i++) {
      let val = matrix[i][lead];
      for (let j = 0; j < colCount; j++) {
        matrix[i][j] -= val * matrix[r][j];
      }
    }
    lead++;
  }
  return matrix;
}

function rref(matrix) {
  if (typeof math.rref === 'function') {
    return math.rref(matrix);
  }
  let m = math.matrix(matrix);
  let nRows = m.size()[0];
  let nCols = m.size()[1];
  let lead = 0;
  let arr = m.toArray();
  for (let r = 0; r < nRows; r++) {
    if (nCols <= lead) return arr;
    let i = r;
    while (arr[i][lead] === 0) {
      i++;
      if (i === nRows) {
        i = r;
        lead++;
        if (nCols === lead) return arr;
      }
    }
    let tmp = arr[i];
    arr[i] = arr[r];
    arr[r] = tmp;

    let val = arr[r][lead];
    if (val !== 0) {
      for (let j = 0; j < nCols; j++) arr[r][j] /= val;
    }

    for (let i = 0; i < nRows; i++) {
      if (i !== r) {
        let val = arr[i][lead];
        for (let j = 0; j < nCols; j++) {
          arr[i][j] -= val * arr[r][j];
        }
      }
    }
    lead++;
  }
  return arr;
}

function powerIteration(A, num_iter=1000, tol=1e-10) {
  const n = A.length;
  let b = Array(n).fill(1);
  let lambda = 0;

  for (let iter = 0; iter < num_iter; iter++) {
    // Compute b_new = A * b
    let b_new = Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        b_new[i] += A[i][j] * b[j];
    // Normalize
    const norm = Math.sqrt(b_new.reduce((sum, x) => sum + x*x, 0));
    b_new = b_new.map(x => x / norm);

    // Rayleigh quotient for eigenvalue estimate
    let num = 0, den = 0;
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let j = 0; j < n; j++)
        s += A[i][j] * b_new[j];
      num += b_new[i] * s;
      den += b_new[i] * b_new[i];
    }
    let lambda_new = num / den;

    // Check for convergence
    if (Math.abs(lambda_new - lambda) < tol)
      break;
    lambda = lambda_new;
    b = b_new;
  }
  return {eigenvalue: lambda, eigenvector: b};
}

function transpose(M) {
  return M[0].map((_, i) => M.map(row => row[i]));
}

function dot(u, v) {
  let sum = 0;
  for (let i = 0; i < u.length; i++) sum += u[i]*v[i];
  return sum;
}

// Gram-Schmidt QR Decomposition
function qrDecomposition(A) {
  let n = A.length;
  let Q = [];
  let R = Array.from({length: n}, () => Array(n).fill(0));

  // Deep copy of A's columns
  let V = [];
  for (let j = 0; j < n; j++)
    V.push(A.map(row => row[j]));

  for (let i = 0; i < n; i++) {
    let v = V[i].slice();
    for (let j = 0; j < i; j++) {
      R[j][i] = dot(Q[j], V[i]);
      for (let k = 0; k < n; k++) v[k] -= R[j][i] * Q[j][k];
    }
    R[i][i] = Math.sqrt(dot(v, v));
    Q.push(v.map(x => x / (R[i][i] || 1e-12)));
  }
  // Reconstruct Q as columns
  let Qmat = Array.from({length: n}, (_, i) => Q.map(col => col[i]));
  return {Q: Qmat, R: R};
}

function qrAlgorithm(A, maxIter=100, tol=1e-8) {
  let n = A.length;
  let Ak = A.map(row => row.slice());
  for (let iter = 0; iter < maxIter; iter++) {
    let {Q, R} = qrDecomposition(Ak);
    // Multiply R*Q
    let newA = Array.from({length: n}, () => Array(n).fill(0));
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        for (let k = 0; k < n; k++)
          newA[i][j] += R[i][k] * Q[k][j];
    // Check for convergence
    let diff = 0;
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        diff += Math.abs(newA[i][j] - Ak[i][j]);
    if (diff < tol) break;
    Ak = newA;
  }
  // Eigenvalues are on diagonal
  return Ak.map((row, i) => row[i]);
}

function findEigenvector(A, lambda) {
  const n = A.length;
  // Construct (A - lambda*I)
  let M = [];
  for (let i = 0; i < n; i++) {
    M.push(A[i].map((val, j) => val - (i === j ? lambda : 0)));
  }
  // Augment with zero column
  for (let i = 0; i < n; i++) M[i].push(0);

  let rrefM = rref(M);

  // Extract eigenvector (free variable = 1)
  // This code picks last variable free, sets to 1
  let x = Array(n).fill(0);
  x[n-1] = 1;
  for (let i = n-2; i >= 0; i--) {
    let sum = 0;
    for (let j = i+1; j < n; j++) {
      sum += rrefM[i][j] * x[j];
    }
    x[i] = -sum;
  }
  return x;
}

function toSqrtOrNumber(x, tol=1e-8) {
  let sq = Math.sqrt(x);
  if (Math.abs(sq - Math.round(sq)) < tol) {
    return `${Math.round(sq)}&sup2;`;
  }
  if (Math.abs(Math.round(sq * sq) - x) < tol) {
    return `&radic;${Math.round(x)}`;
  }
  return null;
}

function formatQuadraticRoots(a, b, c) {
  // a*lambda^2 + b*lambda + c = 0
  if (Math.abs(a) < 1e-12) return [];
  // λ = (-b ± sqrt(b^2-4ac)) / (2a)
  let D = b*b - 4*a*c;
  if (D < 0) {
    return ["Complex eigenvalues"];
  }
  let sqrtD = Math.sqrt(D);
  // Symbolic if not perfect square:
  if (Math.abs(Math.round(sqrtD)*Math.round(sqrtD) - D) < 1e-8) {
    sqrtD = Math.round(sqrtD);
    return [
      `<span>(-${b} + ${sqrtD !== 1 ? sqrtD : ""})/(${2*a})</span>`,
      `<span>(-${b} - ${sqrtD !== 1 ? sqrtD : ""})/(${2*a})</span>`
    ];
  }
  // Otherwise show as ±sqrt
  return [
    `<span>(${(-b)} + &radic;${D})/(${2*a})</span>`,
    `<span>(${(-b)} - &radic;${D})/(${2*a})</span>`
  ];
}

function calculateMatrix() {
  let isAxB = isAxBMode();
  let output = "";
  if (isAxB) {
    const {A, b} = getMatrix();
    output += "<h3>Input Matrix A</h3>" + matrixToHtml(A);
    output += "<h3>Vector b</h3>" + matrixToHtml(b.map(x => [x]));
    const {type, solution, rref, rref: rrefAug} = solveAxEqualsB(A, b);
    output += "<h3>Augmented Matrix REF</h3>" + matrixToHtml(ref(A.map((row, i) => row.concat([b[i]]))));
    output += "<h3>Augmented Matrix RREF</h3>" + matrixToHtml(rrefAug);

    if (type === "unique") {
      output += "<h3>Solution Vector x</h3>" + matrixToHtml(solution.map(x => [x]));
    } else if (type === "none") {
      output += "<h3>No Solution</h3>";
    } else if (type === "infinite") {
      output += "<h3>Infinite Solutions (Free variables)</h3>";
      // Could print rref or more info if you want
    }
  } else {
    // original code
    const matrix = getMatrix();
    output = "<h3>Input Matrix</h3>" + matrixToHtml(matrix);
    const refM = ref(matrix);
    output += "<h3>Row Echelon Form (REF)</h3>" + matrixToHtml(refM);

    const rrefM = rref(matrix);
    output += "<h3>Reduced Row Echelon Form (RREF)</h3>" + matrixToHtml(rrefM);

    if (matrix.length === matrix[0].length) {
      let eigenvalues = qrAlgorithm(matrix, 200, 1e-10);
      output += "<h3>Eigenvalues (approximate)</h3>" + eigenvalues.map(e =>
        e.toFixed(6)
      ).join(', ');

      output += "<h3>Eigenvectors (approximate)</h3>";
      eigenvalues.forEach((lambda, i) => {
        let vec = findEigenvector(matrix, lambda);
        output += `&lambda;<sub>${i+1}</sub> = ${lambda.toFixed(6)}: [${vec.map(x => x.toFixed(6)).join(', ')}]<br>`;
      });   
    }
  }

  document.getElementById('output').innerHTML = output;
}


window.onload = function() {
  onModeChange();
};
