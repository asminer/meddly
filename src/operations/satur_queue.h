/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MEDDLY_SATUR_QUEUE_H
#define MEDDLY_SATUR_QUEUE_H

#include "../defines.h"
#include <vector>

namespace MEDDLY {
    class output;

    /*
        Queue of indexes to explore.
        Quick addition and removal from the queue,
        and also checking if an index belongs to the queue already.
     */
    class index_queue {
            /**
                A linked list of nodes.
                Element i is NOTINQ if it is not in the queue.
                Otherwise, element i is the index of the next item
                in the queue, or -1 if it is the end of the list.
            */
            std::vector <int> data;

            /// Index (in data) of the front of the queue.
            int head;
            /// Index (in data) of the end of the queue.
            int tail;

            static inline int NULPTR() {
                return -1;
            }
            static inline int NOTINQ() {
                return -2;
            }
        public:
            index_queue() {
                head = tail = NULPTR();
            }

            static inline int NOINDEX() {
                return -1;
            }

            inline void resize(unsigned sz) {
                data.resize(sz, NOTINQ());
            }
            inline void clear() {
                data.assign(data.size(), NOTINQ());
            }
            inline bool isEmpty() const {
                return NULPTR() == head;
            }
            inline int front() const {
                return (head >= 0) ? head : NOINDEX();
            }
            inline void add(int i) {
                if (i >= data.size()) {
                    data.resize(i+1, NOTINQ());
                }
                if (NOTINQ() != data[i]) {
                    // already in queue
                    return;
                }
                if (NULPTR() == head) {
                    // queue is empty
                    head = i;
                } else {
                    // queue is not empty
                    MEDDLY_CHECK_RANGE(0, tail, int(data.size()));
                    data[tail] = i;
                }
                tail = i;
                data[i] = NULPTR();
            }
            inline int remove() {
                MEDDLY_CHECK_RANGE(0, head, int(data.size()));
                int ans = head;
                head = data[head];
                data[ans] = NOTINQ();
                return ans;
            }

            void show(output &s) const;
    };

};

#endif
