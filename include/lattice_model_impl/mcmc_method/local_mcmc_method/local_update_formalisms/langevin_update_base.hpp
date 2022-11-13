#ifndef MAIN_LANGEVIN_UPDATE_BASE_UPDATE_HPP
#define MAIN_LANGEVIN_UPDATE_BASE_UPDATE_HPP


#include <lattice_model_impl/mcmc_method/mcmc_method_base.hpp>


namespace lm_impl {
    namespace mcmc_method {
        template<typename MCMCMethod, typename Model>
        class LangevinUpdateBase : public MCMCMethodBase<LangevinUpdateBase<MCMCMethod, Model>> {
        public:
            explicit LangevinUpdateBase(const json params) :
                MCMCMethodBase<LangevinUpdateBase<MCMCMethod, Model>>(params),
                epsilon_(this->template get_entry<double>("epsilon", 0.01)),
                sqrt2epsilon_(sqrt(2 * epsilon_))
            {
                normal_ = std::normal_distribution<double>(0, 1);
            }

            explicit LangevinUpdateBase(const double epsilon=0.01):
                LangevinUpdateBase(json{{"epsilon", epsilon}})
            {}

            template<typename System>
            void initialize(System &system) {
                model_ptr_ = system.get_mcmc_model().get();
            }

            double get_stepsize() const
            {
                return epsilon_;
            }

            template<typename T>
            T estimate_drift_term(const T site) {
                return model_ptr_->get_drift_term(site);
            }

            template<typename T>
            T estimate_drift_term(const T site, const std::vector<T*> neighbours) {
                return model_ptr_->get_drift_term(site, neighbours);
            }

            template<typename T>
            T operator()(const T site) {
                return mcmc_method().update(site, model_ptr_->get_drift_term(site), epsilon_, sqrt2epsilon_);
            }

            template<typename T>
            T operator()(const T site, const std::vector<T*> neighbours) {
                return mcmc_method().update(site, model_ptr_->get_drift_term(site, neighbours), epsilon_, sqrt2epsilon_);
            }

            template<typename T>
            T operator()(const T site, const double epsilon) {
                T eps_drift_term = model_ptr_->get_drift_term(site);
                return mcmc_method().update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
            }

            template<typename T>
            T operator()(const T site, const std::vector<T*> neighbours, const double epsilon) {
                T eps_drift_term = model_ptr_->get_drift_term(site, neighbours);
                return mcmc_method().update(site, eps_drift_term, epsilon, std::sqrt(2 * epsilon));
            }

        protected:           
            Model* model_ptr_;

            std::normal_distribution<double> normal_;

            double epsilon_;
            double sqrt2epsilon_;

        private:
            MCMCMethod &mcmc_method() {
                return *static_cast<MCMCMethod*>(this);
            }

            const MCMCMethod &mcmc_method() const {
                return *static_cast<const MCMCMethod*>(this);
            }
        };
    }
}

#endif //MAIN_LANGEVIN_UPDATE_BASE_UPDATE_HPP
