#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "terrain-mesher/geometry.hpp"

namespace terrain_mesher::core {

template <typename T>
class ObjPool;

template <typename T>
class pool_ptr {
  public:
    pool_ptr() = default;
    pool_ptr(ObjPool<T>* pool, std::size_t index) : pool_(pool), index_(index) {}

    T& operator*() const { return *get(); }
    T* operator->() const { return get(); }
    T* get() const { return pool_->get_addr(index_); }

    ObjPool<T>* pool() const noexcept { return pool_; }
    std::size_t index() const noexcept { return index_; }

    void clear() noexcept {
        pool_ = nullptr;
        index_ = invalid_index_;
    }

    void recycle() {
        if (pool_ != nullptr) {
            pool_->recycle(*this);
            clear();
        }
    }

    bool is_valid() const noexcept { return pool_ != nullptr && pool_->contains(*this); }
    explicit operator bool() const noexcept { return is_valid(); }

    bool operator==(const pool_ptr& other) const noexcept {
        return pool_ == other.pool_ && index_ == other.index_;
    }

    bool operator!=(const pool_ptr& other) const noexcept { return !(*this == other); }

  private:
    static constexpr std::size_t invalid_index_ = ~static_cast<std::size_t>(0);

    ObjPool<T>* pool_ = nullptr;
    std::size_t index_ = invalid_index_;

    template <typename>
    friend class ObjPool;
};

template <typename T>
class ObjPool {
  private:
    struct private_tag {};

  public:
    static std::shared_ptr<ObjPool> create() {
        return std::make_shared<ObjPool>(private_tag{});
    }

    explicit ObjPool(private_tag) {}

    void reserve(std::size_t capacity) { pool_.reserve(capacity); }

    template <typename... Args>
    pool_ptr<T> spawn(Args&&... args) {
        pool_.emplace_back(std::forward<Args>(args)...);
        return pool_ptr<T>(this, pool_.size() - 1);
    }

    void recycle(const pool_ptr<T>&) {}

    bool contains(const pool_ptr<T>& value) const noexcept {
        return value.pool_ == this && value.index_ < pool_.size();
    }

  private:
    T* get_addr(std::size_t index) const { return const_cast<T*>(&pool_[index]); }

    std::vector<T> pool_;

    friend class pool_ptr<T>;
};

class DelaunayTriangle;
class QuadEdge;

using dt_ptr = pool_ptr<DelaunayTriangle>;
using qe_ptr = pool_ptr<QuadEdge>;

class QuadEdge {
  public:
    QuadEdge() = default;

    void init(qe_ptr self_ptr);
    void init(qe_ptr self_ptr, qe_ptr qprev);
    void recycle_next();

    qe_ptr Onext() const noexcept { return next_; }
    qe_ptr Sym() const { return qnext_->qnext_; }
    qe_ptr Rot() const noexcept { return qnext_; }
    qe_ptr invRot() const noexcept { return qprev_; }

    qe_ptr Oprev() const { return Rot()->Onext()->Rot(); }
    qe_ptr Dnext() const { return Sym()->Onext()->Sym(); }
    qe_ptr Dprev() const { return invRot()->Onext()->invRot(); }
    qe_ptr Lnext() const { return invRot()->Onext()->Rot(); }
    qe_ptr Lprev() const { return Onext()->Sym(); }
    qe_ptr Rnext() const { return Rot()->Onext()->invRot(); }
    qe_ptr Rprev() const { return Sym()->Onext(); }

    Point2 Org() const noexcept { return data_; }
    Point2 Dest() const { return Sym()->data_; }

    dt_ptr Lface() const noexcept { return lface_; }
    void set_Lface(dt_ptr value) noexcept { lface_ = value; }

    void set_end_points(const Point2& origin, const Point2& destination);

    friend void splice(qe_ptr a, qe_ptr b);

  private:
    qe_ptr qnext_;
    qe_ptr qprev_;
    qe_ptr self_ptr_;
    Point2 data_;
    qe_ptr next_;
    dt_ptr lface_;
};

void splice(qe_ptr a, qe_ptr b);

class DelaunayTriangle {
  public:
    void init(dt_ptr self_ptr, qe_ptr edge);
    dt_ptr link_to(dt_ptr triangle);
    dt_ptr get_link() const noexcept { return next_face_; }
    qe_ptr get_anchor() const noexcept { return anchor_; }
    void dont_anchor(qe_ptr edge);
    void reshape(qe_ptr edge);

    Point2 point1() const { return anchor_->Org(); }
    Point2 point2() const { return anchor_->Dest(); }
    Point2 point3() const { return anchor_->Lprev()->Org(); }

  private:
    qe_ptr anchor_;
    dt_ptr next_face_;
    dt_ptr self_ptr_;
};

class DelaunayMesh {
  public:
    DelaunayMesh();
    virtual ~DelaunayMesh() = default;

    void init_mesh(const Point2& a, const Point2& b, const Point2& c, const Point2& d);

    virtual bool should_swap(const Point2& point, qe_ptr edge);
    virtual void scan_triangle(dt_ptr triangle);

    qe_ptr locate(const Point2& point);
    qe_ptr locate(const Point2& point, qe_ptr hint);
    qe_ptr spoke(const Point2& point, qe_ptr edge);
    void optimize(const Point2& point, qe_ptr spoke_edge);
    void insert(const Point2& point, dt_ptr triangle);

  protected:
    bool is_interior(qe_ptr edge);

    std::shared_ptr<ObjPool<QuadEdge>> edges_;
    std::shared_ptr<ObjPool<DelaunayTriangle>> triangles_;
    qe_ptr starting_edge_;
    dt_ptr first_face_;

  private:
    dt_ptr make_face(qe_ptr edge);
    void delete_edge(qe_ptr edge);
    qe_ptr connect(qe_ptr a, qe_ptr b);
    void swap(qe_ptr edge);
    bool ccw_boundary(qe_ptr edge);
    bool on_edge(const Point2& point, qe_ptr edge);
    unsigned int next_random_number();

    std::mt19937 random_gen_;
};

inline void QuadEdge::init(qe_ptr self_ptr, qe_ptr qprev) {
    self_ptr_ = self_ptr;
    qprev_ = qprev;
    qprev->qnext_ = self_ptr;
}

inline void QuadEdge::init(qe_ptr self_ptr) {
    self_ptr_ = self_ptr;

    qe_ptr e0 = self_ptr;
    qe_ptr e1 = self_ptr.pool()->spawn();
    qe_ptr e2 = self_ptr.pool()->spawn();
    qe_ptr e3 = self_ptr.pool()->spawn();

    e1->init(e1, e0);
    e2->init(e2, e1);
    e3->init(e3, e2);

    e0->qprev_ = e3;
    e3->qnext_ = e0;

    e0->next_ = e0;
    e1->next_ = e3;
    e2->next_ = e2;
    e3->next_ = e1;
}

inline void QuadEdge::recycle_next() {
    if (!qnext_) {
        return;
    }

    qe_ptr e1 = qnext_;
    qe_ptr e2 = qnext_->qnext_;
    qe_ptr e3 = qprev_;

    e1->qnext_.clear();
    e2->qnext_.clear();
    e3->qnext_.clear();

    e1->recycle_next();
    e1.recycle();
    e2->recycle_next();
    e2.recycle();
    e3->recycle_next();
    e3.recycle();
}

inline void QuadEdge::set_end_points(const Point2& origin, const Point2& destination) {
    data_ = origin;
    Sym()->data_ = destination;
}

inline void splice(qe_ptr a, qe_ptr b) {
    qe_ptr alpha = a->Onext()->Rot();
    qe_ptr beta = b->Onext()->Rot();

    qe_ptr t1 = b->Onext();
    qe_ptr t2 = a->Onext();
    qe_ptr t3 = beta->Onext();
    qe_ptr t4 = alpha->Onext();

    a->next_ = t1;
    b->next_ = t2;
    alpha->next_ = t3;
    beta->next_ = t4;
}

inline void DelaunayTriangle::init(dt_ptr self_ptr, qe_ptr edge) {
    self_ptr_ = self_ptr;
    reshape(edge);
}

inline dt_ptr DelaunayTriangle::link_to(dt_ptr triangle) {
    next_face_ = triangle;
    return self_ptr_;
}

inline void DelaunayTriangle::dont_anchor(qe_ptr edge) {
    if (anchor_ == edge) {
        anchor_ = edge->Lnext();
    }
}

inline void DelaunayTriangle::reshape(qe_ptr edge) {
    anchor_ = edge;
    edge->set_Lface(self_ptr_);
    edge->Lnext()->set_Lface(self_ptr_);
    edge->Lprev()->set_Lface(self_ptr_);
}

inline DelaunayMesh::DelaunayMesh()
    : edges_(ObjPool<QuadEdge>::create()),
      triangles_(ObjPool<DelaunayTriangle>::create()),
      random_gen_(42) {
    edges_->reserve(4096);
    triangles_->reserve(1024);
}

inline dt_ptr DelaunayMesh::make_face(qe_ptr edge) {
    dt_ptr triangle = triangles_->spawn();
    triangle->init(triangle, edge);
    first_face_ = triangle->link_to(first_face_);
    return triangle;
}

inline void DelaunayMesh::init_mesh(const Point2& a,
                                    const Point2& b,
                                    const Point2& c,
                                    const Point2& d) {
    qe_ptr ea = edges_->spawn();
    ea->init(ea);
    ea->set_end_points(a, b);

    qe_ptr eb = edges_->spawn();
    eb->init(eb);
    splice(ea->Sym(), eb);
    eb->set_end_points(b, c);

    qe_ptr ec = edges_->spawn();
    ec->init(ec);
    splice(eb->Sym(), ec);
    ec->set_end_points(c, d);

    qe_ptr ed = edges_->spawn();
    ed->init(ed);
    splice(ec->Sym(), ed);
    ed->set_end_points(d, a);
    splice(ed->Sym(), ea);

    qe_ptr diag = edges_->spawn();
    diag->init(diag);
    splice(ed->Sym(), diag);
    splice(eb->Sym(), diag->Sym());
    diag->set_end_points(a, c);

    starting_edge_ = ea;
    first_face_.clear();

    make_face(ea->Sym());
    make_face(ec->Sym());
}

inline void DelaunayMesh::delete_edge(qe_ptr edge) {
    splice(edge, edge->Oprev());
    splice(edge->Sym(), edge->Sym()->Oprev());
    edge->recycle_next();
    edge.recycle();
}

inline qe_ptr DelaunayMesh::connect(qe_ptr a, qe_ptr b) {
    qe_ptr edge = edges_->spawn();
    edge->init(edge);

    splice(edge, a->Lnext());
    splice(edge->Sym(), b);
    edge->set_end_points(a->Dest(), b->Org());
    return edge;
}

inline void DelaunayMesh::swap(qe_ptr edge) {
    dt_ptr f1 = edge->Lface();
    dt_ptr f2 = edge->Sym()->Lface();

    qe_ptr a = edge->Oprev();
    qe_ptr b = edge->Sym()->Oprev();

    splice(edge, a);
    splice(edge->Sym(), b);
    splice(edge, a->Lnext());
    splice(edge->Sym(), b->Lnext());
    edge->set_end_points(a->Dest(), b->Dest());

    f1->reshape(edge);
    f2->reshape(edge->Sym());
}

inline bool DelaunayMesh::ccw_boundary(qe_ptr edge) {
    return !right_of(edge->Oprev()->Dest(), edge->Org(), edge->Dest());
}

inline bool DelaunayMesh::on_edge(const Point2& point, qe_ptr edge) {
    const double t1 = (point - edge->Org()).length();
    const double t2 = (point - edge->Dest()).length();
    if (t1 < kEpsilon || t2 < kEpsilon) {
        return true;
    }

    const double t3 = (edge->Org() - edge->Dest()).length();
    if (t1 > t3 || t2 > t3) {
        return false;
    }

    Line line(edge->Org(), edge->Dest());
    return std::fabs(line.eval(point)) < kEpsilon;
}

inline bool DelaunayMesh::should_swap(const Point2& point, qe_ptr edge) {
    qe_ptr t = edge->Oprev();
    return in_circle(edge->Org(), t->Dest(), edge->Dest(), point);
}

inline void DelaunayMesh::scan_triangle(dt_ptr) {}

inline bool DelaunayMesh::is_interior(qe_ptr edge) {
    return edge->Lnext()->Lnext()->Lnext() == edge && edge->Rnext()->Rnext()->Rnext() == edge;
}

inline unsigned int DelaunayMesh::next_random_number() {
    return random_gen_() % std::numeric_limits<unsigned int>::max();
}

inline qe_ptr DelaunayMesh::locate(const Point2& point) {
    return locate(point, starting_edge_);
}

inline qe_ptr DelaunayMesh::locate(const Point2& point, qe_ptr hint) {
    qe_ptr edge = hint;
    double t = tri_area(point, edge->Dest(), edge->Org());

    if (t > 0.0) {
        t = -t;
        edge = edge->Sym();
    }

    while (true) {
        qe_ptr edge_origin = edge->Onext();
        qe_ptr edge_dest = edge->Dprev();

        double to = tri_area(point, edge_origin->Dest(), edge_origin->Org());
        double td = tri_area(point, edge_dest->Dest(), edge_dest->Org());

        if (td > 0.0) {
            if (to > 0.0 || (to == 0.0 && t == 0.0)) {
                starting_edge_ = edge;
                return edge;
            }
            t = to;
            edge = edge_origin;
        } else {
            if (to > 0.0) {
                if (td == 0.0 && t == 0.0) {
                    starting_edge_ = edge;
                    return edge;
                }
                t = td;
                edge = edge_dest;
            } else {
                if (t == 0.0 && !left_of(edge_origin->Dest(), edge->Org(), edge->Dest())) {
                    edge = edge->Sym();
                } else if ((next_random_number() & 1U) == 0U) {
                    t = to;
                    edge = edge_origin;
                } else {
                    t = td;
                    edge = edge_dest;
                }
            }
        }
    }
}

inline qe_ptr DelaunayMesh::spoke(const Point2& point, qe_ptr edge) {
    std::array<dt_ptr, 4> new_faces{};
    int face_index = 0;

    qe_ptr boundary_edge;

    dt_ptr left_face = edge->Lface();
    left_face->dont_anchor(edge);
    new_faces[face_index++] = left_face;

    if (on_edge(point, edge)) {
        if (ccw_boundary(edge)) {
            boundary_edge = edge;
        } else {
            dt_ptr sym_face = edge->Sym()->Lface();
            new_faces[face_index++] = sym_face;
            sym_face->dont_anchor(edge->Sym());

            edge = edge->Oprev();
            delete_edge(edge->Onext());
        }
    }

    qe_ptr base = edges_->spawn();
    base->init(base);
    base->set_end_points(edge->Org(), point);
    splice(base, edge);

    starting_edge_ = base;
    do {
        base = connect(edge, base->Sym());
        edge = base->Oprev();
    } while (edge->Lnext() != starting_edge_);

    if (boundary_edge) {
        delete_edge(boundary_edge);
    }

    base = boundary_edge ? starting_edge_->Rprev() : starting_edge_->Sym();
    do {
        if (face_index > 0) {
            new_faces[--face_index]->reshape(base);
        } else {
            make_face(base);
        }
        base = base->Onext();
    } while (base != starting_edge_->Sym());

    return starting_edge_;
}

inline void DelaunayMesh::optimize(const Point2& point, qe_ptr spoke_edge) {
    qe_ptr start_spoke = spoke_edge;
    qe_ptr spoke = spoke_edge;

    do {
        qe_ptr edge = spoke->Lnext();
        if (is_interior(edge) && should_swap(point, edge)) {
            swap(edge);
        } else {
            spoke = spoke->Onext();
            if (spoke == start_spoke) {
                break;
            }
        }
    } while (true);

    spoke = start_spoke;
    do {
        qe_ptr edge = spoke->Lnext();
        dt_ptr triangle = edge->Lface();
        if (triangle) {
            scan_triangle(triangle);
        }
        spoke = spoke->Onext();
    } while (spoke != start_spoke);
}

inline void DelaunayMesh::insert(const Point2& point, dt_ptr triangle) {
    qe_ptr edge = triangle ? locate(point, triangle->get_anchor()) : locate(point);

    if (point == edge->Org() || point == edge->Dest()) {
        optimize(point, edge);
        return;
    }

    qe_ptr start_spoke = spoke(point, edge);
    if (start_spoke) {
        optimize(point, start_spoke->Sym());
    }
}

}  // namespace terrain_mesher::core
