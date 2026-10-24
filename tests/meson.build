
bug_00 = executable(
    'bug_00',
    ['bug_00.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('bug_00', bug_00)

######################################################################

bug_01 = executable(
    'bug_01',
    'bug_01.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('bug_01', bug_01)

######################################################################

bug_02 = executable(
    'bug_02',
    'bug_02.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('bug_02', bug_02)

######################################################################

bug_03 = executable(
    'bug_03',
    'bug_03.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('bug_03', bug_03)

######################################################################

bug_04 = executable(
    'bug_04',
    'bug_04.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('bug_04', bug_04)

######################################################################

check_xA = executable(
    'check_xA',
    'check_xA.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('check_xA', check_xA)

######################################################################

chk_copy = executable(
    'chk_copy',
    'chk_copy.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_copy', chk_copy)

######################################################################

chk_cross = executable(
    'chk_cross',
    ['chk_cross.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_cross', chk_cross)

######################################################################

chk_ct = executable(
    'chk_ct',
    'chk_ct.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_ct', chk_ct)

######################################################################

chk_iters = executable(
    'chk_iters',
    ['chk_iters.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_iters', chk_iters)

######################################################################

chk_minterms = executable(
    'chk_minterms',
    ['chk_minterms.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_minterms --set --fully', chk_minterms, args: ['--set', '--fully'])
test('chk_minterms --set --quasi', chk_minterms, args: ['--set', '--quasi'])
test('chk_minterms --rel --fully', chk_minterms, args: ['--rel', '--fully'])
test('chk_minterms --rel --quasi', chk_minterms, args: ['--rel', '--quasi'])
test('chk_minterms --rel --ident', chk_minterms, args: ['--rel', '--ident'])

######################################################################

chk_mtcoll = executable(
    'chk_mtcoll',
    ['chk_mtcoll.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_mtcoll --set --fully', chk_mtcoll, args: ['--set', '--fully'])
test('chk_mtcoll --set --quasi', chk_mtcoll, args: ['--set', '--quasi'])
test('chk_mtcoll --rel --fully', chk_mtcoll, args: ['--rel', '--fully'])
test('chk_mtcoll --rel --quasi', chk_mtcoll, args: ['--rel', '--quasi'])
test('chk_mtcoll --rel --ident', chk_mtcoll, args: ['--rel', '--ident'])

######################################################################

chk_reorder = executable(
    'chk_reorder',
    'chk_reorder.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('chk_reorder', chk_reorder)

######################################################################

ddedge = executable(
    'ddedge',
    'ddedge.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('ddedge', ddedge)

######################################################################

domain_list = executable(
    'domain_list',
    'domain_list.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('domain_list', domain_list)

######################################################################

kan_batch = executable(
    'kan_batch',
    ['kan_batch.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kan_batch', kan_batch)

######################################################################

kan_chkrs = executable(
    'kan_chkrs',
    ['kan_chkrs.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kan_chkrs -f --mtb --front', kan_chkrs,
    args: ['-f', '--mtb', '--front'])
test('kan_chkrs -f --mtb --trad',  kan_chkrs,
    args: ['-f', '--mtb', '--trad'])
test('kan_chkrs -f --mtb --sat1',  kan_chkrs,
    args: ['-f', '--mtb', '--sat1'])
test('kan_chkrs -f --mti --trad',  kan_chkrs,
    args: ['-f', '--mti', '--trad'])
test('kan_chkrs -f --mti --sat1',  kan_chkrs,
    args: ['-f', '--mti', '--sat1'])
test('kan_chkrs -f --ev  --trad',  kan_chkrs,
    args: ['-f', '--ev', '--trad'])
test('kan_chkrs -f --ev  --sat1',  kan_chkrs,
    args: ['-f', '--ev', '--sat1'])

test('kan_chkrs -b --mtb --front', kan_chkrs,
    args: ['-b', '--mtb', '--front'])
test('kan_chkrs -b --mtb --trad',  kan_chkrs,
    args: ['-b', '--mtb', '--trad'])
test('kan_chkrs -b --mtb --sat1',  kan_chkrs,
    args: ['-b', '--mtb', '--sat1'])
test('kan_chkrs -b --mti --trad',  kan_chkrs,
    args: ['-b', '--mti', '--trad'])
test('kan_chkrs -b --mti --sat1',  kan_chkrs,
    args: ['-b', '--mti', '--sat1'])
test('kan_chkrs -b --ev  --trad',  kan_chkrs,
    args: ['-b', '--ev', '--trad'])
test('kan_chkrs -b --ev  --sat1',  kan_chkrs,
    args: ['-b', '--ev', '--sat1'])

######################################################################

kan_index = executable(
    'kan_index',
    ['kan_index.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kan_index', kan_index)

######################################################################

kan_io = executable(
    'kan_io',
    ['kan_io.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kan_io', kan_io)

######################################################################

kan_relprod = executable(
    'kan_relprod',
    ['kan_relprod.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kan_relprod', kan_relprod)

######################################################################

kan_show = executable(
    'kan_show',
    ['kan_show.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kan_show', kan_show)

######################################################################

kanban = executable(
    'kanban',
    ['kanban.cc', 'simple_model.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('kanban', kanban)

######################################################################

nodestor = executable(
    'nodestor',
    'nodestor.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('nodestor', nodestor)

######################################################################

nqueens = executable(
    'nqueens',
    'nqueens.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('nqueens', nqueens)

######################################################################

ops_comp = executable(
    'ops_comp',
    ['ops_comp.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_comp --set --MT --int', ops_comp,
    args: ['--set', '--MT', '--int'])
test('ops_comp --set --MT --float', ops_comp,
    args: ['--set', '--MT', '--float'])
test('ops_comp --set --EVp --long', ops_comp,
    args: ['--set', '--EVp', '--long'])

test('ops_comp --rel --MT --int', ops_comp,
    args: ['--rel', '--MT', '--int'], timeout: 60)
test('ops_comp --rel --MT --float', ops_comp,
    args: ['--rel', '--MT', '--float'], timeout: 60)
test('ops_comp --rel --EVp --long', ops_comp,
    args: ['--rel', '--EVp', '--long'], timeout: 60)
test('ops_comp --rel --EVs --float', ops_comp,
    args: ['--rel', '--EVs', '--float'], timeout: 60)

######################################################################

ops_distmin = executable(
    'ops_distmin',
    ['ops_distmin.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_distmin --set --int', ops_distmin, args: ['--set', '--int'])
test('ops_distmin --set --float', ops_distmin, args: ['--set', '--float'])
test('ops_distmin --rel --int', ops_distmin, args: ['--rel', '--int'])
test('ops_distmin --rel --float', ops_distmin, args: ['--rel', '--float'])

######################################################################

ops_divmod = executable(
    'ops_divmod',
    ['ops_divmod.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_divmod --set --MT --int', ops_divmod,
    args: ['--set', '--MT', '--int'])
test('ops_divmod --set --MT --float', ops_divmod,
    args: ['--set', '--MT', '--float'])
test('ops_divmod --set --EVp --long', ops_divmod,
    args: ['--set', '--EVp', '--long'])

test('ops_divmod --rel --MT --int', ops_divmod,
    args: ['--rel', '--MT', '--int'])
test('ops_divmod --rel --MT --float', ops_divmod,
    args: ['--rel', '--MT', '--float'])
test('ops_divmod --rel --EVp --long', ops_divmod,
    args: ['--rel', '--EVp', '--long'])
test('ops_divmod --rel --EVs --float', ops_divmod,
    args: ['--rel', '--EVs', '--float'])

######################################################################

ops_maxmin = executable(
    'ops_maxmin',
    ['ops_maxmin.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_maxmin --set --MT --int', ops_maxmin,
    args: ['--set', '--MT', '--int'])
test('ops_maxmin --set --MT --float', ops_maxmin,
    args: ['--set', '--MT', '--float'])
test('ops_maxmin --set --EVp --long', ops_maxmin,
    args: ['--set', '--EVp', '--long'])

test('ops_maxmin --rel --MT --int', ops_maxmin,
    args: ['--rel', '--MT', '--int'])
test('ops_maxmin --rel --MT --float', ops_maxmin,
    args: ['--rel', '--MT', '--float'])
test('ops_maxmin --rel --EVp --long', ops_maxmin,
    args: ['--rel', '--EVp', '--long'])
test('ops_maxmin --rel --EVs --float', ops_maxmin,
    args: ['--rel', '--EVs', '--float'])

######################################################################

ops_mult = executable(
    'ops_mult',
    ['ops_mult.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_mult --set --MT --int', ops_mult,
    args: ['--set', '--MT', '--int'])
test('ops_mult --set --MT --float', ops_mult,
    args: ['--set', '--MT', '--float'])
test('ops_mult --set --EVp --long', ops_mult,
    args: ['--set', '--EVp', '--long'])

test('ops_mult --rel --MT --int', ops_mult,
    args: ['--rel', '--MT', '--int'])
test('ops_mult --rel --MT --float', ops_mult,
    args: ['--rel', '--MT', '--float'])
test('ops_mult --rel --EVp --long', ops_mult,
    args: ['--rel', '--EVp', '--long'])
test('ops_mult --rel --EVs --float', ops_mult,
    args: ['--rel', '--EVs', '--float'])

######################################################################

ops_plusminus = executable(
    'ops_plusminus',
    ['ops_plusminus.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_plusminus --set --MT --int', ops_plusminus,
    args: ['--set', '--MT', '--int'])
test('ops_plusminus --set --MT --float', ops_plusminus,
    args: ['--set', '--MT', '--float'])
test('ops_plusminus --set --EVp --long', ops_plusminus,
    args: ['--set', '--EVp', '--long'])

test('ops_plusminus --rel --MT --int', ops_plusminus,
    args: ['--rel', '--MT', '--int'])
test('ops_plusminus --rel --MT --float', ops_plusminus,
    args: ['--rel', '--MT', '--float'])
test('ops_plusminus --rel --EVp --long', ops_plusminus,
    args: ['--rel', '--EVp', '--long'])
test('ops_plusminus --rel --EVs --float', ops_plusminus,
    args: ['--rel', '--EVs', '--float'])

######################################################################

ops_prepost = executable(
    'ops_prepost',
    ['ops_prepost.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_prepost --set --MTb --inF --relI', ops_prepost,
    args: ['--set', '--MTb', '--inF', '--relI'])
test('ops_prepost --set --MTb --inF --relQ', ops_prepost,
    args: ['--set', '--MTb', '--inF', '--relQ'])
test('ops_prepost --set --MTb --inF --relF', ops_prepost,
    args: ['--set', '--MTb', '--inF', '--relF'])
test('ops_prepost --set --MTb --inQ --relI', ops_prepost,
    args: ['--set', '--MTb', '--inQ', '--relI'])
test('ops_prepost --set --MTb --inQ --relQ', ops_prepost,
    args: ['--set', '--MTb', '--inQ', '--relQ'])
test('ops_prepost --set --MTb --inQ --relF', ops_prepost,
    args: ['--set', '--MTb', '--inQ', '--relF'])

test('ops_prepost --set --MTi --inF --relI', ops_prepost,
    args: ['--set', '--MTi', '--inF', '--relI'])
test('ops_prepost --set --MTi --inF --relQ', ops_prepost,
    args: ['--set', '--MTi', '--inF', '--relQ'])
test('ops_prepost --set --MTi --inF --relF', ops_prepost,
    args: ['--set', '--MTi', '--inF', '--relF'])

test('ops_prepost --set --EVp --inF --relI', ops_prepost,
    args: ['--set', '--EVp', '--inF', '--relI'])
test('ops_prepost --set --EVp --inF --relQ', ops_prepost,
    args: ['--set', '--EVp', '--inF', '--relQ'])
test('ops_prepost --set --EVp --inF --relF', ops_prepost,
    args: ['--set', '--EVp', '--inF', '--relF'])
test('ops_prepost --set --EVp --inQ --relI', ops_prepost,
    args: ['--set', '--EVp', '--inQ', '--relI'])
test('ops_prepost --set --EVp --inQ --relQ', ops_prepost,
    args: ['--set', '--EVp', '--inQ', '--relQ'])
test('ops_prepost --set --EVp --inQ --relF', ops_prepost,
    args: ['--set', '--EVp', '--inQ', '--relF'])

######################################################################

ops_set = executable(
    'ops_set',
    ['ops_set.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_set --set', ops_set, args: ['--set'])
test('ops_rel --rel', ops_set, args: ['--rel'])

######################################################################

ops_user_un = executable(
    'ops_user_un',
    ['ops_user_un.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_user_un --set --MT --int', ops_user_un,
    args: ['--set', '--MT', '--int'])
test('ops_user_un --set --MT --float', ops_user_un,
    args: ['--set', '--MT', '--float'])
test('ops_user_un --set --EVp --long', ops_user_un,
    args: ['--set', '--EVp', '--long'])

test('ops_user_un --rel --MT --int', ops_user_un,
    args: ['--rel', '--MT', '--int'])
test('ops_user_un --rel --MT --float', ops_user_un,
    args: ['--rel', '--MT', '--float'])
test('ops_user_un --rel --EVp --long', ops_user_un,
    args: ['--rel', '--EVp', '--long'])
test('ops_user_un --rel --EVs --float', ops_user_un,
    args: ['--rel', '--EVs', '--float'])

######################################################################

ops_vmmult = executable(
    'ops_vmmult',
    ['ops_vmmult.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('ops_vmmult --int --vF --mF', ops_vmmult, args: ['--int', '--vF', '--mF'])
test('ops_vmmult --int --vF --mI', ops_vmmult, args: ['--int', '--vF', '--mI'])
test('ops_vmmult --int --vF --mQ', ops_vmmult, args: ['--int', '--vF', '--mQ'])
test('ops_vmmult --int --vQ --mF', ops_vmmult, args: ['--int', '--vQ', '--mF'])
test('ops_vmmult --int --vQ --mI', ops_vmmult, args: ['--int', '--vQ', '--mI'])
test('ops_vmmult --int --vQ --mQ', ops_vmmult, args: ['--int', '--vQ', '--mQ'])

test('ops_vmmult --real --vF --mF', ops_vmmult,
    args: ['--real', '--vF', '--mF'])
test('ops_vmmult --real --vF --mI', ops_vmmult,
    args: ['--real', '--vF', '--mI'])
test('ops_vmmult --real --vF --mQ', ops_vmmult,
    args: ['--real', '--vF', '--mQ'])
test('ops_vmmult --real --vQ --mF', ops_vmmult,
    args: ['--real', '--vQ', '--mF'])
test('ops_vmmult --real --vQ --mI', ops_vmmult,
    args: ['--real', '--vQ', '--mI'])
test('ops_vmmult --real --vQ --mQ', ops_vmmult,
    args: ['--real', '--vQ', '--mQ'])

######################################################################

randtest = executable(
    'randtest',
    ['randtest.cc', 'randomize.cc'],
    dependencies: [ gmp_dep, meddly_dep]
)

test('randtest', randtest)

######################################################################

refcount = executable(
    'refcount',
    'refcount.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('refcount', refcount)

######################################################################

sat_test = executable(
    'sat_test',
    'sat_test.cc',
    dependencies: [ gmp_dep, meddly_dep]
)

test('sat_test', sat_test)

