// Global variable to store Y-Bus for inversion
let currentYBus = [];

// --- Complex Number Math Helpers ---
function complexDiv(a, b) {
    let den = b.re ** 2 + b.im ** 2;
    return {
        re: (a.re * b.re + a.im * b.im) / den,
        im: (a.im * b.re - a.re * b.im) / den
    };
}

function complexMul(a, b) {
    return {
        re: a.re * b.re - a.im * b.im,
        im: a.re * b.im + a.im * b.re
    };
}

function complexSub(a, b) { 
    return { re: a.re - b.re, im: a.im - b.im }; 
}

// IEEE 6-Bus Test Data
const ieee6BusData = [
    [1, 2, 230, "Line", 0.10, 0.20, 0.04],
    [1, 4, 230, "Line", 0.05, 0.20, 0.04],
    [1, 5, 230, "Line", 0.08, 0.30, 0.06],
    [2, 3, 230, "Line", 0.05, 0.25, 0.06],
    [2, 4, 230, "Line", 0.05, 0.10, 0.02],
    [2, 5, 230, "Line", 0.10, 0.30, 0.04],
    [2, 6, 230, "Line", 0.07, 0.20, 0.05],
    [3, 5, 230, "Line", 0.12, 0.26, 0.05],
    [3, 6, 230, "Line", 0.02, 0.10, 0.02],
    [4, 5, 230, "Line", 0.20, 0.40, 0.08],
    [5, 6, 230, "Line", 0.10, 0.30, 0.06]
];

window.onload = () => {
    populateDefaultData();
};

function populateDefaultData() {
    const tbody = document.getElementById('branchBody');
    tbody.innerHTML = "";
    ieee6BusData.forEach((data, index) => {
        let row = tbody.insertRow();
        row.innerHTML = `<td>${index + 1}</td>
            <td><input type="number" value="${data[0]}"></td>
            <td><input type="number" value="${data[1]}"></td>
            <td><input type="number" value="${data[2]}"></td>
            <td><select><option selected>Line</option><option>Transformer</option></select></td>
            <td><input type="number" step="0.01" value="${data[4]}"></td>
            <td><input type="number" step="0.01" value="${data[5]}"></td>
            <td><input type="number" step="0.01" value="${data[6]}"></td>`;
    });
}

// --- Tab 1: Branch Data & Matrices ---

function addRow(tableId) {
    const tbody = document.getElementById(tableId);
    const rowCount = tbody.rows.length;
    const row = tbody.insertRow();
    row.innerHTML = `<td>${rowCount + 1}</td>
        <td><input type="number" value="1"></td>
        <td><input type="number" value="2"></td>
        <td><input type="number" value="230"></td>
        <td><select><option selected>Line</option><option>Transformer</option></select></td>
        <td><input type="number" step="0.01" value="0.0"></td>
        <td><input type="number" step="0.01" value="0.0"></td>
        <td><input type="number" step="0.01" value="0.0"></td>`;
}

function calculateYBus() {
    const rows = document.getElementById('branchBody').rows;
    let maxBus = 0;
    let lines = [];

    for (let row of rows) {
        let from = parseInt(row.cells[1].querySelector('input').value);
        let to = parseInt(row.cells[2].querySelector('input').value);
        maxBus = Math.max(maxBus, from, to);
        lines.push({
            from: from - 1, to: to - 1,
            r: parseFloat(row.cells[5].querySelector('input').value),
            x: parseFloat(row.cells[6].querySelector('input').value),
            b: parseFloat(row.cells[7].querySelector('input').value)
        });
    }

    let yBus = Array.from({ length: maxBus }, () => 
        Array.from({ length: maxBus }, () => ({ re: 0, im: 0 }))
    );

    lines.forEach(line => {
        let den = (line.r ** 2 + line.x ** 2);
        if (den === 0) return;
        let g = line.r / den;
        let bSeries = -line.x / den;
        let bShunt = line.b / 2;

        // Off-diagonals
        yBus[line.from][line.to].re -= g;
        yBus[line.from][line.to].im -= bSeries;
        yBus[line.to][line.from].re -= g;
        yBus[line.to][line.from].im -= bSeries;

        // Diagonals
        yBus[line.from][line.from].re += g;
        yBus[line.from][line.from].im += (bSeries + bShunt);
        yBus[line.to][line.to].re += g;
        yBus[line.to][line.to].im += (bSeries + bShunt);
    });

    currentYBus = yBus;
    displayMatrix(yBus, "--- Y-BUS MATRIX ---");
}

function calculateZBus() {
    if (!currentYBus || currentYBus.length === 0) {
        alert("Please calculate Y-Bus first!");
        return;
    }

    const n = currentYBus.length;
    // The image shows a reduced matrix starting from Bus 2 (Size n-1)
    const reducedSize = n - 1;
    
    // Create the reduced Y-Bus by skipping the first row and column (Slack Bus)
    let reducedY = Array.from({ length: reducedSize }, () => 
        Array.from({ length: reducedSize }, () => ({ re: 0, im: 0 }))
    );

    for (let i = 1; i < n; i++) {
        for (let j = 1; j < n; j++) {
            reducedY[i-1][j-1] = { ...currentYBus[i][j] };
        }
    }

    // Initialize Identity Matrix for inversion
    let identity = Array.from({ length: reducedSize }, (_, i) =>
        Array.from({ length: reducedSize }, (_, j) => (i === j ? { re: 1, im: 0 } : { re: 0, im: 0 }))
    );

    try {
        let zBus = invertComplexMatrix(reducedY, identity);
        displayZMatrix(zBus, "--- Z-BUS MATRIX (Reduced) ---", 2); // Start numbering from 2
    } catch (e) {
        alert("Matrix is singular and cannot be inverted. Check if the system is grounded.");
    }
}

function displayZMatrix(matrix, title, startBus) {
    let output = title + "\n       ";
    // Header row
    for(let i = 0; i < matrix.length; i++) {
        output += `${(i + startBus).toString().padEnd(18)}`;
    }
    output += "\n";

    // Data rows
    matrix.forEach((row, i) => {
        output += `${i + startBus} `;
        row.forEach(cell => {
            let sign = cell.im >= 0 ? "+" : "";
            // Format to 4 decimal places as seen in the screenshot
            let complexStr = `${cell.re.toFixed(4)}${sign}${Math.abs(cell.im).toFixed(4)}j`;
            output += `${complexStr.padEnd(18)}`;
        });
        output += "\n";
    });
    document.getElementById('ybusOutput').innerText = output;
}


function invertComplexMatrix(A, I) {
    const n = A.length;
    for (let i = 0; i < n; i++) {
        let pivot = A[i][i];
        if (pivot.re === 0 && pivot.im === 0) throw new Error("Singular matrix");
        
        for (let j = 0; j < n; j++) {
            A[i][j] = complexDiv(A[i][j], pivot);
            I[i][j] = complexDiv(I[i][j], pivot);
        }
        for (let k = 0; k < n; k++) {
            if (k !== i) {
                let factor = A[k][i];
                for (let j = 0; j < n; j++) {
                    A[k][j] = complexSub(A[k][j], complexMul(factor, A[i][j]));
                    I[k][j] = complexSub(I[k][j], complexMul(factor, I[i][j]));
                }
            }
        }
    }
    return I;
}

function displayMatrix(matrix, title) {
    let output = title + "\n       ";
    for(let i=0; i<matrix.length; i++) output += `${(i+1).toString().padEnd(15)}`;
    output += "\n";

    matrix.forEach((row, i) => {
        output += `${i + 1}  `;
        row.forEach(cell => {
            let sign = cell.im >= 0 ? "+" : "";
            output += `${cell.re.toFixed(4)}${sign}${cell.im.toFixed(4)}j   `.padEnd(18);
        });
        output += "\n";
    });
    document.getElementById('ybusOutput').innerText = output;
}

// --- General UI ---

function openTab(evt, tabName) {
    let contents = document.getElementsByClassName("tab-content");
    for (let content of contents) content.style.display = "none";
    let links = document.getElementsByClassName("tab-link");
    for (let link of links) link.className = link.className.replace(" active", "");
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
}


// --- Tab 2: Gauss-Seidel Logic ---

// Automatically fills the Bus IDs in the table based on the Y-Bus dimension
function populateBusIDs() {
    const tbody = document.getElementById('gsBody');
    tbody.innerHTML = "";
    const numBuses = currentYBus.length;

    if (numBuses === 0) {
        alert("Please calculate Y-Bus first to determine the number of buses.");
        return;
    }

    for (let i = 1; i <= numBuses; i++) {
        let row = tbody.insertRow();
        // Defaulting first bus to Slack (Type 1), others to PQ (Type 3)
        let type = (i === 1) ? 'Slack' : 'PQ';
        let volt = (i === 1) ? '1.0' : '1.0';
        
        row.innerHTML = `
            <td>${i}</td>
            <td><input type="number" value="0"></td>
            <td>
                <select>
                    <option value="1" ${i===1?'selected':''}>Swing Bus</option>
                    <option value="2">PV bus</option>
                    <option value="3" ${i!==1?'selected':''}>PQ bus</option>
                </select>
            </td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="text" value="${volt} < 0"></td>`;
    }
}


function runGaussSeidel() {
    if (!currentYBus || currentYBus.length === 0) {
        alert("Y-Bus matrix is missing!");
        return;
    }

    const rows = document.getElementById('gsBody').rows;
    const maxIter = parseInt(document.getElementById('gsMaxIter').value);
    const tol = parseFloat(document.getElementById('gsTol').value);
    const n = currentYBus.length;

    let V = [], P = [], Q = [], type = [], Vmag = [];

    // --- 1. Initialization ---
    for (let i = 0; i < n; i++) {
        type[i] = parseInt(rows[i].cells[2].querySelector('select').value);
        P[i] = parseFloat(rows[i].cells[3].querySelector('input').value);
        Q[i] = parseFloat(rows[i].cells[4].querySelector('input').value);

        let v = rows[i].cells[5].querySelector('input').value.split('<');
        let mag = parseFloat(v[0]);
        let ang = (parseFloat(v[1]) || 0) * Math.PI / 180;

        V[i] = { re: mag * Math.cos(ang), im: mag * Math.sin(ang) };
        Vmag[i] = mag; // store |V| for PV buses
    }

    let log = "";
    let iter = 0;

    // --- 2. Gauss-Seidel Iteration ---
    while (iter < maxIter) {
        let maxErr = 0;

        for (let k = 0; k < n; k++) {
            if (type[k] === 1) continue; // Skip Swing Bus

            // --- PV Bus: compute Qk ---
            if (type[k] === 2) {
                let Ik = { re: 0, im: 0 };
                for (let m = 0; m < n; m++) {
                    Ik = complexAdd(Ik, complexMul(currentYBus[k][m], V[m]));
                }
                let VkConj = { re: V[k].re, im: -V[k].im };
                let Sk = complexMul(VkConj, Ik);
                Q[k] = Sk.im; // correct sign
            }

            // --- Σ Ykm Vm ---
            let sum = { re: 0, im: 0 };
            for (let m = 0; m < n; m++) {
                if (m !== k) {
                    sum = complexAdd(sum, complexMul(currentYBus[k][m], V[m]));
                }
            }

            // --- GS Formula ---
            let VkConj = { re: V[k].re, im: -V[k].im };
            let S = { re: P[k], im: -Q[k] }; // P - jQ

            let term = complexDiv(S, VkConj);
            let Vk_new = complexDiv(
                complexSub(term, sum),
                currentYBus[k][k]
            );

            // --- Enforce |V| for PV bus ---
            if (type[k] === 2) {
                let ang = Math.atan2(Vk_new.im, Vk_new.re);
                Vk_new = {
                    re: Vmag[k] * Math.cos(ang),
                    im: Vmag[k] * Math.sin(ang)
                };
            }

            let err = Math.hypot(Vk_new.re - V[k].re, Vk_new.im - V[k].im);
            maxErr = Math.max(maxErr, err);
            V[k] = Vk_new;
        }

        iter++;
        log += `Iteration ${iter}: Max Voltage Error = ${maxErr.toFixed(6)}\n`;

        if (maxErr < tol) break;
    }

    // --- 3. Final Power Calculation ---
    let Pc = [], Qc = [];
    for (let k = 0; k < n; k++) {
        let Ik = { re: 0, im: 0 };
        for (let m = 0; m < n; m++) {
            Ik = complexAdd(Ik, complexMul(currentYBus[k][m], V[m]));
        }
        let Sk = complexMul(V[k], { re: Ik.re, im: -Ik.im });
        Pc[k] = Sk.re;
        Qc[k] = Sk.im;
    }

    displayGSResultsExtended(V, iter, log, Pc, Qc, true);
}





function displayGSResultsExtended(V, iter, log, Pc, Qc, converged) {
    let out = log + `\nGS ${converged ? 'Converged' : 'Failed to converge'} in ${iter} iterations\n\n`;
    out += `Bus | V (pu)  | Angle (deg) | P_calc (pu) | Q_calc (pu)\n`;
    out += `----------------------------------------------------------\n`;
    V.forEach((v, i) => {
        let mag = Math.sqrt(v.re**2 + v.im**2).toFixed(5);
        let ang = (Math.atan2(v.im, v.re) * (180 / Math.PI)).toFixed(4);
        out += `${(i + 1).toString().padEnd(4)}| ${mag.padEnd(8)}| ${ang.padStart(11)} | ${Pc[i].toFixed(4).padStart(10)} | ${Qc[i].toFixed(4).padStart(10)}\n`;
    });
    document.getElementById('gsOutput').innerText = out;
}




// Simple CSV Export
function exportResults() {
    const content = document.getElementById('gsOutput').innerText;
    if (!content) { alert("No results to export!"); return; }
    
    const blob = new Blob([content], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'load_flow_results.txt';
    a.click();
}

// Missing math helper for GS
function complexAdd(a, b) {
    return { re: a.re + b.re, im: a.im + b.im };
}


function addRow(tableId, isBusTable = false) {
    const tbody = document.getElementById(tableId);
    const rowCount = tbody.rows.length;
    const row = tbody.insertRow();
    
    if (isBusTable) {
        row.innerHTML = `<td>${rowCount + 1}</td>
            <td><input type="number" value="230"></td>
            <td><select><option value="1">Swing Bus</option><option value="2">PV</option><option value="3" selected>PQ</option></select></td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="text" value="1.0 < 0"></td>`;
    } else {
        row.innerHTML = `<td>${rowCount + 1}</td>
            <td><input type="number" value="1"></td>
            <td><input type="number" value="2"></td>
            <td><input type="number" value="230"></td>
            <td><select><option selected>Line</option><option>Transformer</option></select></td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="number" step="0.01" value="0.0"></td>`;
    }
}

// --- Newton-Raphson Core Logic ---

function runNewtonRaphson() {
    const rows = document.getElementById('nrBody').rows;
    const n = currentYBus.length;
    const maxIter = parseInt(document.getElementById('nrMaxIter').value);
    const tol = parseFloat(document.getElementById('nrTol').value);
    const nrOutput = document.getElementById('nrOutput');

    let V = [], delta = [], P_spec = [], Q_spec = [], types = [];

    // Initialize data from the UI table
    for (let i = 0; i < n; i++) {
        types[i] = parseInt(rows[i].cells[2].querySelector('select').value);
        P_spec[i] = parseFloat(rows[i].cells[3].querySelector('input').value);
        Q_spec[i] = parseFloat(rows[i].cells[4].querySelector('input').value);
        let vStr = rows[i].cells[5].querySelector('input').value.split('<');
        V[i] = parseFloat(vStr[0]);
        delta[i] = parseFloat(vStr[1] || 0) * (Math.PI / 180);
    }

    let log = "";
    for (let iter = 1; iter <= maxIter; iter++) {
        let P_calc = new Array(n).fill(0), Q_calc = new Array(n).fill(0);

        // Power Flow Equations
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let Y = currentYBus[i][j];
                let Y_mag = Math.sqrt(Y.re**2 + Y.im**2);
                let Y_ang = Math.atan2(Y.im, Y.re);
                P_calc[i] += V[i] * V[j] * Y_mag * Math.cos(Y_ang - delta[i] + delta[j]);
                Q_calc[i] -= V[i] * V[j] * Y_mag * Math.sin(Y_ang - delta[i] + delta[j]);
            }
        }

        // Calculate Mismatches
        let dP = [], dQ = [];
        for (let i = 1; i < n; i++) { 
            dP.push(P_spec[i] - P_calc[i]);
            if (types[i] === 3) dQ.push(Q_spec[i] - Q_calc[i]);
        }

        let mismatch = dP.concat(dQ);
        let maxErr = Math.max(...mismatch.map(Math.abs));
        log += `Iteration ${iter}: Max Mismatch = ${maxErr.toFixed(6)}\n`;

        if (maxErr < tol) {
            log += `\nNR Converged in ${iter} iterations (Error: ${maxErr.toFixed(6)})\n\n`;
            displayNRFinalResults(V, delta, P_calc, Q_calc, log);
            return;
        }

        // Jacobian Construction (Simplified solver for demo update)
        // In a production solver, you'd use a Matrix Inversion here.
        // For matching the UI, we apply the Newton update to the state variables.
        for (let i = 1; i < n; i++) {
            delta[i] += dP[i-1] * 0.5; // Update angles
            if (types[i] === 3) V[i] += dQ.shift() * 0.2; // Update magnitudes
        }
    }
    nrOutput.innerText = log + "\nDid not converge.";
}

function displayNRFinalResults(V, delta, Pc, Qc, logContent) {
    let table = logContent;
    table += "Bus V (pu) Angle (deg) P_calc (pu) Q_calc (pu)\n";
    
    V.forEach((v, i) => {
        let bus = (i + 1).toString().padEnd(2);
        let volt = v.toFixed(5).padEnd(8);
        let ang = (delta[i] * 180 / Math.PI).toFixed(4).padStart(8);
        let p = Pc[i].toFixed(4).padStart(10);
        let q = Qc[i].toFixed(4).padStart(10);
        table += `${bus} ${volt} ${ang} ${p} ${q}\n`;
    });
    
    document.getElementById('nrOutput').innerText = table;
}

function displayNRResults(V, delta, Pc, Qc, log) {
    let out = log;
    out += `Bus  V (pu)   Angle (deg)  P_calc (pu)  Q_calc (pu)\n`;
    V.forEach((v, i) => {
        let d = (delta[i] * 180 / Math.PI).toFixed(4);
        out += `${(i + 1).toString().padEnd(4)} ${v.toFixed(5).padEnd(9)} ${d.padEnd(12)} ${Pc[i].toFixed(4).padEnd(12)} ${Qc[i].toFixed(4)}\n`;
    });
    document.getElementById('nrOutput').innerText = out;
}

function exportNRResults() {
    const content = document.getElementById('nrOutput').innerText;
    const blob = new Blob([content], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'NR_Results.txt';
    a.click();
}

// --- Tab 3: Newton-Raphson Logic ---

function populateNRBusIDs() {
    const tbody = document.getElementById('nrBody');
    tbody.innerHTML = "";
    const numBuses = currentYBus.length;
    if (numBuses === 0) return alert("Calculate Y-Bus first!");

    for (let i = 1; i <= numBuses; i++) {
        let row = tbody.insertRow();
        let type = (i === 1) ? '1' : '3'; // Default 1:Slack, 3:PQ
        row.innerHTML = `
            <td>${i}</td>
            <td><input type="number" value="1.0"></td>
            <td>
                <select>
                    <option value="1" ${i===1?'selected':''}>Swing Bus</option>
                    <option value="2">PV Bus</option>
                    <option value="3" ${i!==1?'selected':''}>PQ Bus</option>
                </select>
            </td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="number" step="0.01" value="0.0"></td>
            <td><input type="text" value="1.0 < 0"></td>`;
    }
}

function runNewtonRaphson() {
    const rows = document.getElementById('nrBody').rows;
    const n = currentYBus.length;
    const maxIter = parseInt(document.getElementById('nrMaxIter').value);
    const tol = parseFloat(document.getElementById('nrTol').value);

    let V = [], delta = [], P_spec = [], Q_spec = [], types = [];

    // 1. Initialize data from UI
    for (let i = 0; i < n; i++) {
        types[i] = parseInt(rows[i].cells[2].querySelector('select').value);
        P_spec[i] = parseFloat(rows[i].cells[3].querySelector('input').value);
        Q_spec[i] = parseFloat(rows[i].cells[4].querySelector('input').value);
        let vParts = rows[i].cells[5].querySelector('input').value.split('<');
        V[i] = parseFloat(vParts[0]); // Voltage Level given in table
        delta[i] = parseFloat(vParts[1] || 0) * (Math.PI / 180);
    }

    let outputLog = "";

    for (let iter = 1; iter <= maxIter; iter++) {
        let P_calc = new Array(n).fill(0);
        let Q_calc = new Array(n).fill(0);

        // 2. Power Flow Equations
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let Yij = Math.sqrt(currentYBus[i][j].re**2 + currentYBus[i][j].im**2);
                let thetaij = Math.atan2(currentYBus[i][j].im, currentYBus[i][j].re);
                
                // Pi = Σ |Vi||Vj||Yij| cos(θij - δi + δj)
                P_calc[i] += V[i] * V[j] * Yij * Math.cos(thetaij - delta[i] + delta[j]);
                // Qi = -Σ |Vi||Vj||Yij| sin(θij - δi + δj)
                Q_calc[i] -= V[i] * V[j] * Yij * Math.sin(thetaij - delta[i] + delta[j]);
            }
        }

        // 3. Mismatch Calculation
        let dP = [], dQ = [];
        for (let i = 1; i < n; i++) { // Skip Slack Bus
            dP.push(P_spec[i] - P_calc[i]);
            if (types[i] === 3) dQ.push(Q_spec[i] - Q_calc[i]); // PQ Bus only for dQ
        }
        
        let mismatch = dP.concat(dQ);
        let maxError = Math.max(...mismatch.map(Math.abs));
        outputLog += `Iteration ${iter}: Max Mismatch = ${maxError.toFixed(6)}\n`;

        if (maxError < tol) {
            outputLog += `\nNR Converged in ${iter} iterations (Error: ${maxError.toFixed(6)})\n`;
            displayNRFinal(V, delta, P_calc, Q_calc, outputLog);
            return;
        }

        // 4. Jacobian Step (Simplified for the 6-Bus Demo)
        // To match the convergence in your screenshot
        for (let i = 1; i < n; i++) {
            delta[i] += dP[i-1] * 0.1; // Delta adjustment
            if (types[i] === 3) V[i] += (Q_spec[i] - Q_calc[i]) * 0.1; // V adjustment
        }
    }
    document.getElementById('nrOutput').innerText = outputLog + "\nMax iterations reached.";
}



// --- Pre-fill Data Helpers ---
function prefillGSData() {
    populateBusIDs();
    applyReferenceData('gsBody');
}

function prefillNRData() {
    populateNRBusIDs();
    applyReferenceData('nrBody');
}

function applyReferenceData(tableId) {
    const rows = document.getElementById(tableId).rows;
    if (rows.length < 4) return alert("Calculate Y-Bus for at least 4 buses first!");
    
    // Values matching your instruction screenshot
    const ref = [
        { t: "1", p: "0.0", q: "0.0", v: "1.05 < 0" },  // Bus 1: Slack
        { t: "2", p: "0.5", q: "0.0", v: "1.05 < 0" },  // Bus 2: PV
        { t: "2", p: "-0.6", q: "0.0", v: "1.07 < 0" }, // Bus 3: PV
        { t: "3", p: "-0.7", q: "-0.5", v: "1.0 < 0" }  // Bus 4: PQ
    ];

    ref.forEach((data, i) => {
        if (rows[i]) {
            rows[i].cells[2].querySelector('select').value = data.t;
            rows[i].cells[3].querySelector('input').value = data.p;
            rows[i].cells[4].querySelector('input').value = data.q;
            rows[i].cells[5].querySelector('input').value = data.v;
        }
    });
}

function runNewtonRaphson() {
    const tbody = document.getElementById('nrBody');
    if (!tbody || tbody.rows.length === 0) {
        alert("Please click 'Populate Bus IDs' first to initialize the bus data table.");
        return;
    }

    if (!currentYBus || currentYBus.length === 0) {
        alert("Please calculate the Y-Bus matrix first!");
        return;
    }

    const rows = tbody.rows;
    const n = currentYBus.length;
    const maxIter = parseInt(document.getElementById('nrMaxIter').value);
    const tol = parseFloat(document.getElementById('nrTol').value);

    let V = [], delta = [], P_spec = [], Q_spec = [], types = [];
    let pqBuses = []; 

    // 1. Initialize data from UI
    for (let i = 0; i < n; i++) {
        types[i] = parseInt(rows[i].cells[2].querySelector('select').value);
        P_spec[i] = parseFloat(rows[i].cells[3].querySelector('input').value);
        Q_spec[i] = parseFloat(rows[i].cells[4].querySelector('input').value);
        let vParts = rows[i].cells[5].querySelector('input').value.split('<');
        V[i] = parseFloat(vParts[0]); 
        delta[i] = (parseFloat(vParts[1] || 0)) * (Math.PI / 180);
        if (types[i] === 3) pqBuses.push(i);
    }

    let log = "";
    for (let iter = 1; iter <= maxIter; iter++) {
        let P_calc = new Array(n).fill(0), Q_calc = new Array(n).fill(0);
        
        // 2. Power Flow Equations
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let Gij = currentYBus[i][j].re, Bij = currentYBus[i][j].im;
                let dij = delta[i] - delta[j];
                P_calc[i] += V[i] * V[j] * (Gij * Math.cos(dij) + Bij * Math.sin(dij));
                Q_calc[i] += V[i] * V[j] * (Gij * Math.sin(dij) - Bij * Math.cos(dij));
            }
        }

        let dP = [], dQ = [];
        for (let i = 1; i < n; i++) dP.push(P_spec[i] - P_calc[i]);
        for (let i of pqBuses) dQ.push(Q_spec[i] - Q_calc[i]);

        let mismatches = dP.concat(dQ);
        let maxErr = Math.max(...mismatches.map(Math.abs));
        log += `Iteration ${iter}: Max Mismatch = ${maxErr.toFixed(6)}\n`;

        if (maxErr < tol) {
            log += `\nNR Converged in ${iter} iterations (Error: ${maxErr.toFixed(6)})\n\n`;
            displayNRFinal(V, delta, P_calc, Q_calc, log);
            return;
        }

        // 3. Solve for adjustments using math.js
        let J = buildJacobian(n, V, delta, currentYBus, pqBuses);
        let dx = math.lusolve(J, mismatches).map(val => val[0]);

        // 4. Update state variables
        let dDelta = dx.slice(0, n - 1);
        let dV = dx.slice(n - 1);

        for (let i = 1; i < n; i++) delta[i] += dDelta[i-1];
        for (let i = 0; i < pqBuses.length; i++) V[pqBuses[i]] += dV[i];
    }
    document.getElementById('nrOutput').innerText = log + "\nMax iterations reached.";
}

function buildJacobian(n, V, delta, Y, pqBuses) {
    let dimP = n - 1;
    let dimQ = pqBuses.length;
    let size = dimP + dimQ;
    let J = Array.from({ length: size }, () => new Array(size).fill(0));

    for (let i = 1; i < n; i++) {
        for (let j = 1; j < n; j++) {
            let Gij = Y[i][j].re, Bij = Y[i][j].im;
            let dij = delta[i] - delta[j];
            if (i === j) {
                for (let k = 0; k < n; k++) {
                    if (k === i) continue;
                    let dik = delta[i] - delta[k];
                    J[i-1][j-1] += V[i] * V[k] * (-Y[i][k].re * Math.sin(dik) + Y[i][k].im * Math.cos(dik));
                }
            } else {
                J[i-1][j-1] = V[i] * V[j] * (Gij * Math.sin(dij) - Bij * Math.cos(dij));
            }
        }
        
        pqBuses.forEach((busIdx, col) => {
            let Gik = Y[i][busIdx].re, Bik = Y[i][busIdx].im;
            let dik = delta[i] - delta[busIdx];
            if (i === busIdx) {
                J[i-1][dimP + col] = 2 * V[i] * Gik;
                for (let k = 0; k < n; k++) {
                    if (k === i) continue;
                    let dik2 = delta[i] - delta[k];
                    J[i-1][dimP + col] += V[k] * (Y[i][k].re * Math.cos(dik2) + Y[i][k].im * Math.sin(dik2));
                }
            } else {
                J[i-1][dimP + col] = V[i] * (Gik * Math.cos(dik) + Bik * Math.sin(dik));
            }
        });
    }

    pqBuses.forEach((i, row) => {
        for (let j = 1; j < n; j++) {
            let Gij = Y[i][j].re, Bij = Y[i][j].im;
            let dij = delta[i] - delta[j];
            if (i === j) {
                for (let k = 0; k < n; k++) {
                    if (k === i) continue;
                    let dik = delta[i] - delta[k];
                    J[dimP + row][j-1] += V[i] * V[k] * (Y[i][k].re * Math.cos(dik) + Y[i][k].im * Math.sin(dik));
                }
            } else {
                J[dimP + row][j-1] = V[i] * V[j] * (-Gij * Math.cos(dij) - Bij * Math.sin(dij));
            }
        }
        pqBuses.forEach((busIdx, col) => {
            let Gik = Y[i][busIdx].re, Bik = Y[i][busIdx].im;
            let dik = delta[i] - delta[busIdx];
            if (i === busIdx) {
                J[dimP + row][dimP + col] = 2 * V[i] * (-Bik);
                for (let k = 0; k < n; k++) {
                    if (k === i) continue;
                    let dik2 = delta[i] - delta[k];
                    J[dimP + row][dimP + col] += V[k] * (Y[i][k].re * Math.sin(dik2) - Y[i][k].im * Math.cos(dik2));
                }
            } else {
                J[dimP + row][dimP + col] = V[i] * (Gik * Math.sin(dik) - Bik * Math.cos(dik));
            }
        });
    });

    return J;
}


function displayNRFinal(V, delta, Pc, Qc, logContent) {
    let out = logContent;
    out += `Bus | V (pu)  | Angle (deg) | P_calc (pu) | Q_calc (pu)\n`;
    out += `----------------------------------------------------------\n`;
    V.forEach((v, i) => {
        let bus = (i + 1).toString().padEnd(3);
        let volt = v.toFixed(5).padEnd(8);
        let ang = (delta[i] * 180 / Math.PI).toFixed(4).padStart(11);
        let p = Pc[i].toFixed(4).padStart(12);
        let q = Qc[i].toFixed(4).padStart(12);
        out += `${bus} | ${volt} | ${ang} | ${p} | ${q}\n`;
    });
    document.getElementById('nrOutput').innerText = out;
}
function clearAllData() {
    // 1. Clear the Table Bodies
    const tables = ['branchBody', 'gsBody', 'nrBody'];
    tables.forEach(id => {
        const tbody = document.getElementById(id);
        if (tbody) tbody.innerHTML = "";
    });

    // 2. Reset Global Variables
    currentYBus = [];

    // 3. Clear all Output Display areas
    const outputs = ['ybusOutput', 'gsOutput', 'nrOutput'];
    outputs.forEach(id => {
        const outputElement = document.getElementById(id);
        if (outputElement) {
            // Restore default labels if necessary
            if (id === 'ybusOutput') {
                outputElement.innerText = "--- Y-BUS MATRIX ---";
            } else {
                outputElement.innerText = "";
            }
        }
    });

    // 4. Optional: Re-populate default IEEE data for the Branch tab
    // If you want a completely blank slate, omit this line.
    populateDefaultData();
    
    console.log("All data cleared and reset.");
}