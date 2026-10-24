
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

test('kan_chkrs -f --mtb --front', kan_chkrs, args: ['-f', '--mtb', '--front'])
test('kan_chkrs -f --mtb --trad',  kan_chkrs, args: ['-f', '--mtb', '--trad'])
test('kan_chkrs -f --mtb --sat1',  kan_chkrs, args: ['-f', '--mtb', '--sat1'])
test('kan_chkrs -f --mti --trad',  kan_chkrs, args: ['-f', '--mti', '--trad'])
test('kan_chkrs -f --mti --sat1',  kan_chkrs, args: ['-f', '--mti', '--sat1'])
test('kan_chkrs -f --ev  --trad',  kan_chkrs, args: ['-f', '--ev', '--trad'])
test('kan_chkrs -f --ev  --sat1',  kan_chkrs, args: ['-f', '--ev', '--sat1'])

test('kan_chkrs -b --mtb --front', kan_chkrs, args: ['-b', '--mtb', '--front'])
test('kan_chkrs -b --mtb --trad',  kan_chkrs, args: ['-b', '--mtb', '--trad'])
test('kan_chkrs -b --mtb --sat1',  kan_chkrs, args: ['-b', '--mtb', '--sat1'])
test('kan_chkrs -b --mti --trad',  kan_chkrs, args: ['-b', '--mti', '--trad'])
test('kan_chkrs -b --mti --sat1',  kan_chkrs, args: ['-b', '--mti', '--sat1'])
test('kan_chkrs -b --ev  --trad',  kan_chkrs, args: ['-b', '--ev', '--trad'])
test('kan_chkrs -b --ev  --sat1',  kan_chkrs, args: ['-b', '--ev', '--sat1'])

