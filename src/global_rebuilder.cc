// $Id$

#include "defines.h"
#include "meddly.h"
#include "meddly_expert.h"
#include "hash_stream.h"

size_t MEDDLY::global_rebuilder::RestrictKeyHasher::operator()(const RestrictKey &key) const
{
	MEDDLY::hash_stream s;
	s.start(0);
	s.push(key.p, key.var, key.idx);
	return s.finish();
}

size_t MEDDLY::global_rebuilder::TransformKeyHasher::operator()(const TransformKey &key) const
{
	MEDDLY::hash_stream s;
	s.start(0);
	s.push(key.p, key.var);
	return s.finish();
}

MEDDLY::global_rebuilder::global_rebuilder(expert_forest* source, expert_forest* target)
	: _source(source), _target(target)
{
	if(_source->getDomain() != _target->getDomain()) {
		throw error::DOMAIN_MISMATCH;
	}
}

MEDDLY::dd_edge MEDDLY::global_rebuilder::rebuild(const dd_edge& e)
{
	if(e.getForest() != _source) {
		throw error(error::FOREST_MISMATCH);
	}

	int num_var = e.getForest()->getDomain()->getNumVariables();
	node_handle p = e.getNode();
	node_handle pt = transform(p, num_var);

	dd_edge out(_target);
	out.set(pt);
	return out;
}

MEDDLY::node_handle MEDDLY::global_rebuilder::transform(node_handle p, int target_level)
{
	if(_source->isTerminalNode(p)) {
		return p;
	}

	int top_var = _target->getVarByLevel(target_level);
	while(isLevelAbove(_source->getLevelByVar(top_var), _source->getNodeLevel(p))){
		target_level--;
		top_var = _target->getVarByLevel(target_level);
	}

	auto search = _computed_transform.find({p, top_var});
	if(search != _computed_transform.end()) {
		return _target->linkNode(search->second);
	}

	int size = _target->getVariableSize(top_var);
	node_builder& nb = _target->useNodeBuilder(target_level, size);
	for(int i=0; i<size; i++) {
		node_handle pr = restrict(p, top_var, i);
//		_computed_restrict.clear();
		_source->unlinkNode(pr);
		nb.d(i) = transform(pr, target_level-1);
	}
	node_handle pt = _target->createReducedNode(-1, nb);

	// Saved in the computed table
	_computed_transform.insert({{p, top_var}, pt});

//	_source->linkNode(p);
	return pt;
}

MEDDLY::node_handle MEDDLY::global_rebuilder::restrict(node_handle p, int var, int idx)
{
	int level1 = _source->getNodeLevel(p);
	int level2 = _source->getLevelByVar(var);
	if(level1 < level2) {
		return _source->linkNode(p);
	}
	else {
		auto search = _computed_restrict.find({p, var, idx});
		if(search != _computed_restrict.end()) {
			return _source->linkNode(search->second);
		}

		if(level1 > level2) {
			node_reader* nr = _source->initNodeReader(p, true);

			int size = _source->getVariableSize(_source->getVarByLevel(level1));
			node_builder& nb = _source->useNodeBuilder(level1, size);
			for(int i=0; i<size; i++) {
				nb.d(i) = restrict(nr->d(i), var, idx);
			}
			node_reader::recycle(nr);

			node_handle pr = _source->createReducedNode(-1, nb);

			// Saved in the computed table
			_computed_restrict.insert({{p, var, idx}, pr});

			return _source->linkNode(pr);
//			return pr;
		}
		else {
			node_handle d = _source->getDownPtr(p, idx);
			return _source->linkNode(d);
		}
	}
}

void MEDDLY::global_rebuilder::clearCache()
{
	for(auto cache : _computed_restrict) {
		_source->unlinkNode(cache.second);
	}
	_computed_restrict.clear();
	_computed_transform.clear();
}
