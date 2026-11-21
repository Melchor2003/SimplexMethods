/*SOLVER_START*/
const EPS = 1e-9;

function parseExpr(str){
  const s = String(str).replace(/\s+/g,'').replace(/[−–]/g,'-');
  if(!s) throw new Error('An expression is empty.');
  const re = /([+-]?)(\d*\.?\d+|\d+\.?)?\*?([a-zA-Z][a-zA-Z0-9_]*)/y;
  const out = {}; let pos = 0;
  while(pos < s.length){
    re.lastIndex = pos;
    const m = re.exec(s);
    if(!m) throw new Error('Cannot read "' + s.slice(pos) + '". Use terms like 3x1 + 5x2.');
    const c = (m[2] === undefined ? 1 : parseFloat(m[2])) * (m[1] === '-' ? -1 : 1);
    out[m[3]] = (out[m[3]] || 0) + c;
    pos = re.lastIndex;
  }
  return out;
}

function toFrac(x){
  if(Math.abs(x) < 1e-10) return '0';
  const sgn = x < 0 ? '-' : ''; const ax = Math.abs(x);
  let h1=1,h0=0,k1=0,k0=1,b=ax;
  for(let i=0;i<24;i++){
    const a = Math.floor(b);
    [h1,h0] = [a*h1+h0, h1]; [k1,k0] = [a*k1+k0, k1];
    if(Math.abs(ax - h1/k1) < 1e-9) break;
    b = 1/(b-a); if(!isFinite(b)) break;
  }
  if(k1 > 10000 || Math.abs(ax - h1/k1) > 1e-7) return dec(x);
  return sgn + (k1 === 1 ? h1 : h1 + '/' + k1);
}
function dec(x){
  if(Math.abs(x) < 1e-10) return '0';
  return String(parseFloat(x.toFixed(4)));
}

function solveLP({type, objective, constraints}){
  const names = [];
  const add = o => Object.keys(o).forEach(k => { if(!names.includes(k)) names.push(k); });
  add(objective); constraints.forEach(c => add(c.coefs));
  const n = names.length, m = constraints.length;
  const sign = type === 'max' ? 1 : -1;
  const cost = names.map(k => sign * (objective[k] || 0));

  const rows = constraints.map(c => {
    let co = names.map(k => c.coefs[k] || 0), op = c.op, rhs = c.rhs;
    if(rhs < 0){ co = co.map(v => -v); rhs = -rhs; op = op === '<=' ? '>=' : op === '>=' ? '<=' : '='; }
    return {co, op, rhs};
  });

  const cols = names.map(nm => ({name:nm, type:'x'}));
  const extra = [];
  rows.forEach((r,i) => {
    if(r.op === '<=') extra.push({name:'s'+(i+1), type:'slack', row:i, coef:1, basic:true});
    else if(r.op === '>='){
      extra.push({name:'s'+(i+1), type:'surplus', row:i, coef:-1});
      extra.push({name:'a'+(i+1), type:'art', row:i, coef:1, basic:true});
    } else extra.push({name:'a'+(i+1), type:'art', row:i, coef:1, basic:true});
  });
  extra.forEach(e => cols.push({name:e.name, type:e.type}));
  const total = cols.length;

  const A = rows.map(r => { const a = new Array(total+1).fill(0); r.co.forEach((v,j) => a[j] = v); a[total] = r.rhs; return a; });
  const basis = new Array(m).fill(-1);
  extra.forEach((e,idx) => { A[e.row][n+idx] = e.coef; if(e.basic) basis[e.row] = n+idx; });

  const zm = new Array(total+1).fill(0), zc = new Array(total+1).fill(0);
  for(let j=0;j<n;j++) zc[j] = -cost[j];
  cols.forEach((c,j) => { if(c.type === 'art') zm[j] = 1; });
  for(let i=0;i<m;i++){
    if(cols[basis[i]].type === 'art'){ for(let j=0;j<=total;j++) zm[j] -= A[i][j]; }
  }

  const A0 = A.map(r => r.slice());
  const snap = () => ({ A:A.map(r => r.map(v => Math.abs(v)<1e-12?0:v)), zm:zm.map(v => Math.abs(v)<1e-12?0:v), zc:zc.map(v => Math.abs(v)<1e-12?0:v), basis:basis.slice() });
  const isNeg = j => zm[j] < -EPS || (Math.abs(zm[j]) <= EPS && zc[j] < -EPS);
  const less = (a,b) => Math.abs(zm[a]-zm[b]) > EPS ? zm[a] < zm[b] : zc[a] < zc[b];

  const steps = []; let status = 'optimal', iter = 0;
  while(true){
    let k = -1;
    for(let j=0;j<total;j++){
      if(!isNeg(j)) continue;
      if(k === -1) k = j; else if(iter < 50 && less(j,k)) k = j;
    }
    if(k === -1) break;
    if(iter++ >= 200){ status = 'iterlimit'; break; }
    const ratios = A.map(r => r[k] > EPS ? r[total]/r[k] : null);
    let r = -1;
    ratios.forEach((v,i) => { if(v !== null && (r === -1 || v < ratios[r] - 1e-12 || (Math.abs(v-ratios[r]) <= 1e-12 && basis[i] < basis[r]))) r = i; });
    const st = {table:snap(), enter:k, leave:r, ratios};
    steps.push(st);
    if(r === -1){ status = 'unbounded'; break; }
    const p = A[r][k];
    A[r] = A[r].map(v => v/p);
    for(let i=0;i<m;i++){ if(i !== r && A[i][k] !== 0){ const f = A[i][k]; A[i] = A[i].map((v,j) => v - f*A[r][j]); } }
    const fm = zm[k], fc = zc[k];
    for(let j=0;j<=total;j++){ zm[j] -= fm*A[r][j]; zc[j] -= fc*A[r][j]; }
    basis[r] = k;
    st.pivot = p;
  }

  const finalTable = snap();
  if(status === 'optimal'){
    const artBad = basis.some((b,i) => cols[b].type === 'art' && A[i][total] > 1e-7);
    if(artBad || Math.abs(zm[total]) > 1e-7) status = 'infeasible';
  }

  const values = {}; names.forEach(nm => values[nm] = 0);
  basis.forEach((b,i) => { if(b < n) values[names[b]] = Math.abs(A[i][total]) < 1e-12 ? 0 : A[i][total]; });
  const Z = sign * zc[total];

  let multiple = false;
  if(status === 'optimal'){
    for(let j=0;j<total;j++){
      if(cols[j].type !== 'art' && !basis.includes(j) && Math.abs(zm[j]) < EPS && Math.abs(zc[j]) < EPS) multiple = true;
    }
  }

  // standard form text
  const term = (v, label, first) => {
    const a = Math.abs(v), body = (a === 1 ? '' : toFrac(a)) + label;
    return first ? (v < 0 ? '-' : '') + body : (v < 0 ? ' - ' : ' + ') + body;
  };
  const lines = A0.map((row,i) => {
    let s = '', first = true;
    for(let j=0;j<total;j++){ if(Math.abs(row[j]) > EPS){ s += term(row[j], cols[j].name, first); first = false; } }
    return s + ' = ' + toFrac(row[total]);
  });
  let objExpr = '';
  { let first = true; names.forEach(nm => { const v = objective[nm] || 0; if(Math.abs(v) > EPS){ objExpr += term(v, nm, first); first = false; } }); if(first) objExpr = '0'; }
  let objStr = (type === 'max' ? 'Max Z = ' : 'Min Z = ') + objExpr;
  const arts = cols.filter(c => c.type === 'art');
  let zEq = type === 'max' ? 'z' : "z'";
  names.forEach((nm,j) => { const v = -cost[j]; if(Math.abs(v) > EPS) zEq += term(v, nm, false); });
  arts.forEach(a => { zEq += ' + M' + a.name; });
  zEq += ' = 0';
  const stdObj = type === 'max'
    ? objStr + arts.map(a => ' - M' + a.name).join('')
    : objStr + '   →   Max Z\' = -Z' + arts.map(a => ' - M' + a.name).join('');

  return {status, steps, finalTable, cols, names, values, Z, type, multiple, std:{obj:stdObj, lines, objExpr, zEq}, model:{type, objective, constraints}, total, m};
}

/* numbers behind the "Pivotal Elimination Steps" tables */
function elimData(t, enter, leave){
  const p = t.A[leave][enter], inv = 1/p;
  const newPiv = t.A[leave].map(v => v*inv);
  const ops = [{kind:'pivot', row:leave, old:t.A[leave].slice(), factor:inv, res:newPiv}];
  t.A.forEach((row,i) => {
    if(i === leave) return;
    const f = row[enter];
    ops.push({kind:'row', row:i, f, old:row.slice(), res:row.map((v,j) => v - f*newPiv[j])});
  });
  const fm = t.zm[enter], fc = t.zc[enter];
  ops.push({kind:'z', fm, fc, oldm:t.zm.slice(), oldc:t.zc.slice(),
    resm:t.zm.map((v,j) => v - fm*newPiv[j]), resc:t.zc.map((v,j) => v - fc*newPiv[j])});
  return {ops, newPiv};
}
/*SOLVER_END*/

/* ---------- UI ---------- */
const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];
const esc = s => String(s).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const sub = n => esc(n).replace(/^([a-zA-Z]+)(\d+)$/, '$1<sub>$2</sub>');
let fracMode = true, lastRes = null;
const fmt = x => fracMode ? toFrac(x) : dec(x);
function fmtPair(m, c){
  if(Math.abs(m) < 1e-10) return fmt(c);
  const ms = (Math.abs(m) === 1 ? '' : fmt(Math.abs(m))) + 'M';
  if(Math.abs(c) < 1e-10) return (m < 0 ? '-' : '') + ms;
  return fmt(c) + (m > 0 ? ' + ' : ' - ') + ms;
}

/* theme */
(function(){
  let t = null;
  try{ t = localStorage.getItem('sx_theme'); }catch(e){}
  if(!t) t = matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  const apply = v => { document.documentElement.setAttribute('data-theme', v); $('#theme').textContent = v === 'dark' ? '☀️' : '🌙'; try{ localStorage.setItem('sx_theme', v); }catch(e){} };
  apply(t);
  $('#theme').addEventListener('click', () => apply(document.documentElement.getAttribute('data-theme') === 'dark' ? 'light' : 'dark'));
})();

/* routing */
const titles = {home:'Simplex Linear Programming Solver', solver:'Simplex Method Calculator', team:'Our Team', reviewer:'Reviewer & Downloads'};
function route(){
  let h = (location.hash || '#/home').replace('#/','');
  if(!titles[h]) h = 'home';
  Object.keys(titles).forEach(r => {
    $('#page-'+r).hidden = r !== h;
    const a = document.querySelector('[data-route="'+r+'"]');
    a.classList.toggle('active', r === h);
    if(r === h) a.setAttribute('aria-current','page'); else a.removeAttribute('aria-current');
  });
  document.title = titles[h];
  window.scrollTo(0,0);
  $('#links').classList.remove('open'); $('#burger').setAttribute('aria-expanded','false');
}
addEventListener('hashchange', route);
$('#burger').addEventListener('click', () => { const o = $('#links').classList.toggle('open'); $('#burger').setAttribute('aria-expanded', o); });

/* solver form */
const examples = {
  cls:{t:'max', obj:'120x + 100y', c:[['2x + 2y','<=',8],['5x + 3y','<=',15]]},
  max:{t:'max', obj:'3x1 + 5x2', c:[['x1','<=',4],['2x2','<=',12],['3x1 + 2x2','<=',18]]},
  min:{t:'min', obj:'2x1 + 3x2', c:[['x1 + x2','>=',4],['x1 + 3x2','>=',6]]},
  mix:{t:'max', obj:'2x1 + 3x2 + x3', c:[['x1 + x2 + x3','=',10],['x1 - x2','<=',2],['x2 + x3','<=',8]]},
  unb:{t:'max', obj:'x1 + x2', c:[['x1 - x2','<=',2]]},
  inf:{t:'max', obj:'x1 + x2', c:[['x1 + x2','<=',2],['x1 + x2','>=',5]]}
};
function addConstraint(expr='', op='<=', rhs=''){
  const d = document.createElement('div'); d.className = 'crow';
  d.innerHTML = '<input class="in c-expr" placeholder="e.g. 2x1 + x2" autocomplete="off" spellcheck="false" aria-label="Constraint left side">' +
    '<select class="in c-op" aria-label="Operator"><option value="<=">≤</option><option value=">=">≥</option><option value="=">=</option></select>' +
    '<input class="in c-rhs" type="number" step="any" placeholder="RHS" aria-label="Right side value">' +
    '<button class="btn danger sm c-del" type="button">Remove</button>';
  d.querySelector('.c-expr').value = expr; d.querySelector('.c-op').value = op; d.querySelector('.c-rhs').value = rhs;
  d.querySelector('.c-del').addEventListener('click', () => { if($$('.crow').length > 1){ d.remove(); chips(); } else showErr('Keep at least one constraint.'); });
  d.addEventListener('input', chips);
  $('#cons').appendChild(d); chips();
}
function chips(){
  const names = [];
  const grab = v => { try{ Object.keys(parseExpr(v)).forEach(k => names.includes(k) || names.push(k)); }catch(e){} };
  grab($('#obj').value); $$('.c-expr').forEach(i => grab(i.value));
  $('#nonneg').innerHTML = '<span class="nn-lock" title="Always on">🔒 Non-negativity · always active</span>' +
    (names.length ? names.map(n => '<span class="chip">' + sub(n) + ' ≥ 0</span>').join('') : '<span class="chip dim">all variables ≥ 0</span>');
}
function showErr(msg){ const e = $('#err'); e.textContent = msg; e.hidden = !msg; }
const exNames = {cls:'Maximize 120x + 100y (class example)', max:'Maximize 3x1 + 5x2', min:'Minimize 2x1 + 3x2 (≥ constraints)', mix:'Maximize 2x1 + 3x2 + x3 (≤ and = constraints)', unb:'Unbounded example', inf:'Infeasible example'};
function loadExample(k){
  const ex = examples[k]; if(!ex) return;
  $(ex.t === 'max' ? '#o-max' : '#o-min').checked = true; zlabel();
  $('#obj').value = ex.obj; $('#cons').innerHTML = '';
  ex.c.forEach(c => addConstraint(c[0], c[1], c[2]));
  clearBad(); showErr(''); chips();
  $('#exNote').textContent = 'Loaded: ' + exNames[k] + '. Press Solve LP Model.';
  $('#example').value = '';
}
function zlabel(){ $('#zlab').textContent = ($('#o-max').checked ? 'Max' : 'Min') + ' Z ='; }
function resetForm(){ $('#obj').value = ''; $('#cons').innerHTML = ''; addConstraint(); addConstraint(); $('#example').value = ''; $('#exNote').textContent = ''; clearBad(); showErr(''); $('#results').hidden = true; lastRes = null; chips(); }

function markBad(el){ el.classList.add('invalid'); }
function clearBad(){ $$('.invalid').forEach(e => e.classList.remove('invalid')); }
function readModel(){
  clearBad();
  const type = $('#o-max').checked ? 'max' : 'min';
  let objective;
  try{ objective = parseExpr($('#obj').value); }catch(e){ markBad($('#obj')); throw new Error('Objective function: ' + e.message); }
  if(!Object.values(objective).some(v => Math.abs(v) > EPS)){ markBad($('#obj')); throw new Error('The objective function needs at least one non-zero coefficient.'); }
  const constraints = $$('.crow').map((row,i) => {
    const ex = row.querySelector('.c-expr'), rh = row.querySelector('.c-rhs'), raw = ex.value.trim(), rv = rh.value;
    if(!raw){ markBad(ex); throw new Error('Constraint ' + (i+1) + ' is missing its left side.'); }
    let coefs;
    try{ coefs = parseExpr(raw); }catch(e){ markBad(ex); throw new Error('Constraint ' + (i+1) + ': ' + e.message); }
    if(!Object.values(coefs).some(v => Math.abs(v) > EPS)){ markBad(ex); throw new Error('Constraint ' + (i+1) + ' has only zero coefficients.'); }
    if(rv === '' || !isFinite(parseFloat(rv))){ markBad(rh); throw new Error('Constraint ' + (i+1) + ' needs a numeric right-hand side.'); }
    return {coefs, op: row.querySelector('.c-op').value, rhs: parseFloat(rv)};
  });
  const vars = new Set(Object.keys(objective)); constraints.forEach(c => Object.keys(c.coefs).forEach(k => vars.add(k)));
  if(vars.size > 10) throw new Error('Please use at most 10 variables.');
  return {type, objective, constraints};
}
function solveClick(){
  showErr('');
  let res;
  try{
    clearBad();
    if(!$('#obj').value.trim()){ markBad($('#obj')); throw new Error('Please enter an objective function.'); }
    res = solveLP(readModel());
  }catch(e){ showErr(e.message); $('#results').hidden = true; return; }
  lastRes = res; render();
  $('#results').hidden = false;
  $('#results').scrollIntoView({behavior:'smooth', block:'start'});
}

function tableHTML(t, cols, opt){
  const total = cols.length, {enter=-1, leave=-1, ratios=null} = opt || {};
  let h = '<div class="tscroll"><table class="tab"><thead><tr><th>Basis</th>';
  cols.forEach((c,j) => h += '<th class="' + (j === enter ? 'pcol' : '') + '">' + sub(c.name) + '</th>');
  h += '<th>RHS</th>' + (ratios ? '<th>Ratio</th>' : '') + '</tr></thead><tbody>';
  t.A.forEach((row,i) => {
    h += '<tr class="' + (i === leave ? 'prow' : '') + '"><th>' + sub(cols[t.basis[i]].name) + '</th>';
    for(let j=0;j<total;j++) h += '<td class="' + (j === enter ? 'pcol' : '') + (i === leave && j === enter ? ' pivot' : '') + '">' + fmt(row[j]) + '</td>';
    h += '<td>' + fmt(row[total]) + '</td>';
    if(ratios) h += '<td class="' + (i === leave ? 'minr' : '') + '">' + (ratios[i] === null ? '—' : fmt(ratios[i])) + '</td>';
    h += '</tr>';
  });
  h += '<tr class="zrow"><th>Z</th>';
  for(let j=0;j<total;j++) h += '<td class="' + (j === enter ? 'pcol' : '') + '">' + fmtPair(t.zm[j], t.zc[j]) + '</td>';
  h += '<td>' + fmtPair(t.zm[total], t.zc[total]) + '</td>' + (ratios ? '<td></td>' : '') + '</tr></tbody></table></div>';
  return h;
}

function joinTerms(items){
  return items.map((t,i) => {
    const neg = t.c < 0, body = fmt(Math.abs(t.c)) + '(' + fmt(t.v) + ')';
    return i === 0 ? (neg ? '-' : '') + body : (neg ? ' - ' : ' + ') + body;
  }).join('');
}
function elimHTML(r, s){
  if(s.leave === -1) return '';
  const t = s.table, d = elimData(t, s.enter, s.leave), total = r.cols.length;
  const en = sub(r.cols[s.enter].name);
  const heads = d.ops.map(op => {
    if(op.kind === 'pivot') return sub(r.cols[t.basis[op.row]].name) + '(1/P) = ' + en;
    if(op.kind === 'row'){ const nm = sub(r.cols[t.basis[op.row]].name); return nm + ' - ' + en + '(P' + nm + ') = ' + nm; }
    return 'z - ' + en + '(Pz) = z';
  });
  let h = '<div class="elim"><h4>Pivotal Elimination Steps</h4><div class="tscroll"><table class="tab et"><thead><tr><th>Col</th>' + heads.map(x => '<th>' + x + '</th>').join('') + '</tr></thead><tbody>';
  for(let j=0;j<=total;j++){
    h += '<tr><th>' + (j < total ? sub(r.cols[j].name) : 'RHS') + '</th>';
    d.ops.forEach(op => {
      let txt;
      if(op.kind === 'pivot') txt = fmt(op.old[j]) + '(' + fmt(op.factor) + ') = ' + fmt(op.res[j]);
      else if(op.kind === 'row') txt = fmt(op.old[j]) + ' - ' + fmt(d.newPiv[j]) + '(' + fmt(op.f) + ') = ' + fmt(op.res[j]);
      else txt = fmtPair(op.oldm[j], op.oldc[j]) + ' - ' + fmt(d.newPiv[j]) + '(' + fmtPair(op.fm, op.fc) + ') = ' + fmtPair(op.resm[j], op.resc[j]);
      h += '<td>' + txt + '</td>';
    });
    h += '</tr>';
  }
  return h + '</tbody></table></div><p class="hint">First column: 1/P means 1 ÷ pivot element. Other columns: P followed by a row name is that row\'s entry in the pivot column.</p></div>';
}
function checkHTML(r){
  if(r.status !== 'optimal') return '';
  const v = r.values, tol = 1e-7, sym = {'<=':'≤','>=':'≥','=':'='};
  const okb = '<span class="badge ok">✓ satisfied</span>', badb = '<span class="badge bad">✗ violated</span>';
  let rows = r.model.constraints.map((c,i) => {
    const items = r.names.filter(n => Math.abs(c.coefs[n] || 0) > EPS).map(n => ({c:c.coefs[n], v:v[n]}));
    const lhs = items.reduce((a,t) => a + t.c*t.v, 0);
    const pass = c.op === '<=' ? lhs <= c.rhs + tol : c.op === '>=' ? lhs >= c.rhs - tol : Math.abs(lhs - c.rhs) <= tol;
    return '<tr><th>Constraint ' + (i+1) + '</th><td class="l">' + joinTerms(items) + ' = ' + fmt(lhs) + '</td><td>' + sym[c.op] + ' ' + fmt(c.rhs) + '</td><td>' + (pass ? okb : badb) + '</td></tr>';
  }).join('');
  const nn = r.names.every(n => v[n] >= -1e-9);
  rows += '<tr><th>Non-negativity</th><td class="l">' + r.names.map(n => sub(n) + ' = ' + fmt(v[n])).join(', ') + '</td><td>all ≥ 0</td><td>' + (nn ? okb : badb) + '</td></tr>';
  const zItems = r.names.filter(n => Math.abs(r.model.objective[n] || 0) > EPS).map(n => ({c:r.model.objective[n], v:v[n]}));
  const zc = zItems.reduce((a,t) => a + t.c*t.v, 0);
  rows += '<tr><th>Objective Z</th><td class="l">' + joinTerms(zItems) + ' = ' + fmt(zc) + '</td><td>Z = ' + fmt(r.Z) + '</td><td>' + (Math.abs(zc - r.Z) <= 1e-6 ? '<span class="badge ok">✓ matches</span>' : badb) + '</td></tr>';
  return '<div class="panel" style="margin-top:18px"><h3>Solution Check</h3><p class="hint" style="margin:4px 0 12px">The optimal values are substituted back into the original model.</p><div class="tscroll"><table class="tab ck"><thead><tr><th>Check</th><th>Substitution</th><th>Compare</th><th>Result</th></tr></thead><tbody>' + rows + '</tbody></table></div></div>';
}

function render(){
  if(!lastRes) return;
  const r = lastRes, out = $('#out');
  const badge = {optimal:['ok','✅ Optimal solution found'], unbounded:['warn','⚠️ Unbounded problem'], infeasible:['bad','⛔ No feasible solution'], iterlimit:['warn','⚠️ Iteration limit reached']}[r.status];
  let h = '<div class="summary"><div class="panel"><span class="badge ' + badge[0] + '">' + badge[1] + '</span>';
  if(r.status === 'optimal'){
    h += '<div class="zval">' + (r.type === 'max' ? 'Max' : 'Min') + ' Z = ' + fmt(r.Z) + '</div>' +
      '<small style="color:var(--muted)">after ' + r.steps.length + ' iteration' + (r.steps.length === 1 ? '' : 's') + '</small><div class="vars">' +
      r.names.map(n => '<div class="var"><span>' + sub(n) + '</span><b>' + fmt(r.values[n]) + '</b></div>').join('') + '</div>';
    if(r.multiple) h += '<div class="note">Alternative optima exist: a non-basic variable has a zero Z-row value in the final tableau.</div>';
  } else if(r.status === 'unbounded'){
    h += '<p style="margin-top:14px;color:var(--muted)">The entering variable has no positive entries in its column, so the objective can improve forever. Check that your constraints bound the problem.</p>';
  } else if(r.status === 'infeasible'){
    h += '<p style="margin-top:14px;color:var(--muted)">An artificial variable stays in the final basis with a positive value, so the constraints contradict each other.</p>';
  } else h += '<p style="margin-top:14px;color:var(--muted)">The solver stopped after 200 iterations. The model may be cycling.</p>';
  h += '</div><div class="panel"><h3>Standard form</h3>' +
    '<div class="lpm"><div class="lpm-t">LP Model</div><div class="lpm-row"><span class="lpm-a">' + esc(r.std.objExpr) + '</span><span class="lpm-arrow">--&gt;</span><span class="lpm-z">' + esc(r.std.zEq) + '</span></div></div>' +
    '<div class="std">' + esc(r.std.obj) + '\n' + r.std.lines.map(esc).join('\n') + '</div>' +
    '<div class="nn"><span class="nn-lock">🔒 Non-negativity · always active</span><div class="chips">' + r.cols.map(c => '<span class="chip">' + sub(c.name) + ' ≥ 0</span>').join('') + '</div></div>' +
    (r.cols.some(c => c.type === 'art') ? '<p class="hint">Artificial variables (a) are penalized with the Big M method. M is a very large positive number.</p>' : '') + '</div></div>';

  r.steps.forEach((s,i) => {
    h += '<details class="step" open><summary><span class="num">' + (i+1) + '</span>Tableau ' + (i+1) + (i === 0 ? ' · initial' : '') + '<span class="chev">▾</span></summary><div class="body">' +
      tableHTML(s.table, r.cols, {enter:s.enter, leave:s.leave, ratios:s.ratios});
    h += '<div class="expl"><div><b>Entering variable</b>' + sub(r.cols[s.enter].name) + ' · most negative Z-row value (' + fmtPair(s.table.zm[s.enter], s.table.zc[s.enter]) + ')</div>';
    if(s.leave === -1) h += '<div><b>Leaving variable</b>None: every entry in the pivot column is ≤ 0, so the problem is unbounded.</div>';
    else h += '<div><b>Leaving variable</b>' + sub(r.cols[s.table.basis[s.leave]].name) + ' · smallest ratio (' + fmt(s.ratios[s.leave]) + ')</div><div><b>Pivot element</b>' + fmt(s.table.A[s.leave][s.enter]) + '</div>';
    h += '</div>' + elimHTML(r, s) + '</div></details>';
  });
  if(r.status !== 'unbounded'){
    h += '<details class="step" open><summary><span class="num">' + (r.steps.length + 1) + '</span>' + 'Tableau ' + (r.steps.length + 1) + (r.status === 'optimal' ? ' · final (optimal)' : ' · final') + '<span class="chev">▾</span></summary><div class="body">' + tableHTML(r.finalTable, r.cols) + '</div></details>';
  }
  h += checkHTML(r);
  h += '<div class="legend"><span><i style="background:var(--pivot)"></i>Pivot element</span><span><i style="background:var(--pcol)"></i>Entering column</span><span><i style="background:var(--prow)"></i>Leaving row</span></div>';
  if(r.type === 'min') h += '<p class="hint" style="margin-top:10px">Minimization is solved as Max Z\' = −Z, so the Z row above tracks Z\'. The reported Z is converted back.</p>';
  out.innerHTML = h;
}

/* team */

const IMG_FOLDER = 'images/';
const IMG_TYPES = ['jpg','jpeg','png','webp','gif'];
const PROF = {name:'PTL. Norcelito Galvan', photo:'prof'};
const team = [
  {n:'Rupha Alenna O. Lopega',    r:'Documentation',      l:'Gamay, Northern Samar',    q:'Coding is the art of turning ideas into reality.', photo:'rupha'},
  {n:'Shiena M. Lucban',          r:'UI/UX Designer',     l:'UEP Main Campus',          q:'Design creates culture.',                          photo:'Shiena'},
  {n:'Jethro G. Lucindo',         r:'Quality Assurance',  l:'Gamay,Northern Samar',     q:'Quality is never an accident.',                    photo:'jethro'},
  {n:'Keith Ablow B. Mabanag',    r:'Back-End Developer', l:'UEP Main Campus',          q:'Words give clarity to ideas.',                     photo:'keith'},
  {n:'Miguel Angelo R. Mercado',  r:'Mathematician',      l:'Catarman',                 q:'Numbers never lie.',                               photo:'miguel'},
  {n:'Andrea Paula C. Mopas',     r:'Mathematician',      l:'UEP Main Campus',          q:'One Step at a time.',                              photo:'paula'},
  {n:'Melchor M. Morada Jr.',     r:'UI/UX Designer',     l:'Las Navas',                q:'If you die, I will kill you!.',                    photo:'melchor'},
  {n:'Aljhen L. Nabus',           r:'Documentation',      l:'Pambujan',                 q:'When the time is right I the Lord make it happen.',photo:'aljhen'}
];
function avatar(el, photo, name){
  const initials = name.replace(/^PTL\.\s*/,'').split(/\s+/).map(w => w[0]).slice(0,2).join('').toUpperCase();
  const fallback = () => { el.innerHTML = '<div class="ph">' + esc(initials) + '</div>'; };
  if(!photo){ fallback(); return; }
  /* a photo value that already has an extension (rocky.jpg) is used as written */
  const hasExt = /\.[a-z0-9]{3,4}$/i.test(photo);
  const tries = hasExt ? [photo] : IMG_TYPES.map(t => photo + '.' + t);
  let k = 0;
  const img = new Image(); img.alt = name; img.loading = 'lazy';
  img.onerror = () => { k++; if(k < tries.length) img.src = IMG_FOLDER + tries[k]; else fallback(); };
  img.src = IMG_FOLDER + tries[0]; el.appendChild(img);
}
function buildTeam(){
  const g = $('#teamGrid');
  team.forEach(m => {
    const c = document.createElement('div'); c.className = 'card hover member';
    c.innerHTML = '<div class="avatar"></div><h3>' + esc(m.n) + '</h3><span class="role">' + esc(m.r) + '</span><div class="loc">📍 ' + esc(m.l) + '</div><p class="quote">“' + esc(m.q) + '”</p>';
    avatar(c.querySelector('.avatar'), m.photo, m.n); g.appendChild(c);
  });
  avatar($('#profAvatar'), PROF.photo, PROF.name);
  avatar($('#profAvatar2'), PROF.photo, PROF.name);
}

/* modules (downloads) */
/* ============================================================
   HOW TO ADD YOUR MODULE FILES
   1. Make a folder named  assets  next to index.html
   2. Put your PDF / PPT files in it
   3. Write each file's exact name (with extension) in  file  below
   To add a module, copy one line. To remove one, delete its line.
   ============================================================ */
const FILE_FOLDER = 'assets/';
const MODULES = [
  {tag:'Module 1', title:'Linear Programming', file:'Module-1-Quantitative-Methods.pdf',
   desc:'Intro to Linear Programming: inequalities, graphical method, feasible region, decision variables, objective functions, and formulating LP models.'},
  {tag:'Module 2', title:'Simplex Method', file:'Module-2-QUANTITATIVE-METHODS2.pdf',
   desc:'Covers simplex algorithm steps, slack/surplus/artificial variables, tableaus, pivot operations, and solving LP problems beyond graphing.'},
  {tag:'Module 3', title:'Two-Phase Method & M-Method', file:'Module-3-M-METHOD_AND_TP_Final.pdf',
   desc:'Covers the complete procedure for solving Linear Programming Problems using the Big M Method and the Two-Phase Method, including rules, steps, tableaux, and detailed examples.'},
  {tag:'Module 4', title:'Dual Method', file:'Module-4-The-Dual-Method.pdf',
   desc:'Introduction to Duality in Linear Programming: dual formulation, primal–dual relationships, and solving LP using dual method principles.'},
  {tag:'Module 5', title:'Additional Materials', file:'qm-LP_090138.pptx',
   desc:'Supplementary reviewer and lecture slide for Linear Programming & Simplex topics.'}
];
function buildModules(){
  const g = $('#modulesGrid'); if(!g) return;
  g.innerHTML = MODULES.map(m => {
    const ext = (m.file.split('.').pop() || '').toUpperCase();
    const icon = /^(PPT|PPTX)$/.test(ext) ? '📊' : '📄';
    return '<div class="card hover module"><span class="tag">' + esc(m.tag) + '</span><h3>📘 ' + esc(m.title) + '</h3><p>' + esc(m.desc) + '</p>' +
      '<a class="btn primary sm" href="' + esc(FILE_FOLDER + m.file) + '" download target="_blank" rel="noopener">' + icon + ' Download ' + esc(m.tag) + ' (' + esc(ext) + ')</a></div>';
  }).join('');
}

/* feedback */
let rating = 0, mem = [];
const ratingTxt = ['Select a rating','Poor','Fair','Good','Very Good','Excellent'];
function loadFb(){ try{ return JSON.parse(localStorage.getItem('sx_feedback') || '[]'); }catch(e){ return mem; } }
function saveFb(a){ mem = a; try{ localStorage.setItem('sx_feedback', JSON.stringify(a)); }catch(e){} }
function stars(){
  const box = $('#stars'); box.innerHTML = '';
  for(let i=1;i<=5;i++){
    const b = document.createElement('button'); b.type = 'button'; b.textContent = '★'; b.className = i <= rating ? 'on' : '';
    b.setAttribute('role','radio'); b.setAttribute('aria-checked', i === rating); b.setAttribute('aria-label', i + ' star' + (i>1?'s':'') + ' - ' + ratingTxt[i]);
    b.addEventListener('click', () => { rating = i; stars(); }); box.appendChild(b);
  }
  $('#starTxt').textContent = ratingTxt[rating];
}
function renderFb(){
  const a = loadFb();
  $('#fbList').innerHTML = a.length ? a.slice().reverse().map(f =>
    '<div class="fb"><div class="top"><span>' + esc(f.name) + '</span><span class="st">' + '★'.repeat(f.rating) + '</span></div><p>' + esc(f.text) + '</p><small>' + esc(f.date) + '</small></div>').join('')
    : '<p class="hint" style="margin-top:10px">No feedback yet. Be the first to share your thoughts.</p>';
}
$('#fbForm').addEventListener('submit', e => {
  e.preventDefault();
  const name = $('#fbName').value.trim(), text = $('#fbText').value.trim(), er = $('#fbErr');
  if(!name || !rating || !text){ er.textContent = 'Please add your name, a rating, and your feedback.'; er.hidden = false; return; }
  er.hidden = true;
  const a = loadFb(); a.push({name, rating, text, date:new Date().toLocaleDateString(undefined,{year:'numeric',month:'short',day:'numeric'})});
  saveFb(a); $('#fbName').value = ''; $('#fbText').value = ''; rating = 0; stars(); renderFb();
});

/* wire up */
$('#solve').addEventListener('click', solveClick);
$('#reset').addEventListener('click', resetForm);
$('#addc').addEventListener('click', () => addConstraint());
$('#obj').addEventListener('input', chips);
document.addEventListener('input', e => { if(e.target.classList) e.target.classList.remove('invalid'); });
$('#obj').addEventListener('keydown', e => { if(e.key === 'Enter') solveClick(); });
$$('input[name="otype"]').forEach(i => i.addEventListener('change', zlabel));
$('#example').addEventListener('change', e => loadExample(e.target.value));
$('#f-frac').addEventListener('click', () => { fracMode = true; $('#f-frac').classList.add('on'); $('#f-dec').classList.remove('on'); render(); });
$('#f-dec').addEventListener('click', () => { fracMode = false; $('#f-dec').classList.add('on'); $('#f-frac').classList.remove('on'); render(); });

buildTeam(); buildModules(); stars(); renderFb(); addConstraint(); addConstraint(); route();
