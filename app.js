function generateMatrix() {
  const m = parseInt(document.getElementById('rows').value);
  const n = parseInt(document.getElementById('cols').value);
  const form = document.getElementById('matrix-form');
  form.innerHTML = '';
  form.style.gridTemplateColumns = `repeat(${n}, 70px)`;
  for (let i = 0; i < m; ++i) {
    for (let j = 0; j < n; ++j) {
      const inp = document.createElement('input');
      inp.type = 'number';
      inp.step = 'any';
      inp.id = `cell-${i}-${j}`;
      inp.value = (i === j) ? 1 : 0;
      inp.inputMode = "numeric";
      inp.pattern = "\\d*";
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
  for (let i = 0; i < m; ++i) {
    let row = [];
    for (let j = 0; j < n; ++j) {
      let val = parseFloat(document.getElementById(`cell-${i}-${j}`).value);
      row.push(isNaN(val) ? 0 : val);
    }
    matrix.push(row);
  }
  return matrix;
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


window.onload = generateMatrix;
