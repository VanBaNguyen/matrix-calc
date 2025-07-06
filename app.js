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

  // Dynamically set the bracket height/content
  setMatrixBrackets(m);
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

function calculateMatrix() {
  const matrix = getMatrix();
  let output = "<h3>Input Matrix</h3>" + matrixToHtml(matrix);

  const refM = ref(matrix);
  output += "<h3>Row Echelon Form (REF)</h3>" + matrixToHtml(refM);

  const rrefM = rref(matrix);
  output += "<h3>Reduced Row Echelon Form (RREF)</h3>" + matrixToHtml(rrefM);

  document.getElementById('output').innerHTML = output;
}

window.onload = generateMatrix;
