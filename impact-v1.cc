/* */


#include <impact.h>

#include <aspect/utilities.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <cmath>
#include <iostream>

namespace aspect
{
  namespace InitialTemperature
  {
    //no dim declare
    template <int dim>
    double
    Impact<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // calculate the isobaric core radius
      const double iso_core_radius =    radius_impactor * std::pow(10, -0.346) * std::pow(velocity_impactor, 0.211);
      const double density = reference_density_value;
      const double specific_heat = reference_specific_heat_value;
      //first, get the Coordinate information(spherical or box) and calculate separately
      //second, calculate the distance between all points and impact position
      Point<dim> impact_point;
      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
        {
          const double outer_radius = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()).outer_radius();
          if (dim==2)
            {
              impact_point(0) = outer_radius,
              impact_point(1) = impact_position_one; //question: whether the distance function is good in spherical?
            }
          else if (dim==3)
            {
              impact_point(0) = outer_radius,
              impact_point(1) = impact_position_one,
              impact_point(2) = impact_position_two;
            }
        }

      double distance_to_impact = std::abs(impact_point.distance(position));
      //third, look up material properties and calculate the impact temperature for each points
      // typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, this->n_compositional_fields());
      // typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, this->n_compositional_fields());
      if (distance_to_impact<iso_core_radius)
        {
          double Pshock = density * (intercept_EOS*1000 + slope_EOS*particle_velocity*1000)*particle_velocity*1000;
          double beta = density * std::pow(intercept_EOS*1000,2.0)/(2*slope_EOS);
          double f = -( Pshock / beta ) *std::pow(( 1 - std::sqrt( 2*Pshock/beta + 1) ), -1);
          double temperature_impact = std::pow(specific_heat,-1) * ((Pshock/(2*density*slope_EOS))  *(1-std::pow(f,-1)) - (f-std::log(f)-1)*std::pow(intercept_EOS*1000/slope_EOS, 2.0)); 
        }
        else if (distance_to_impact>=iso_core_radius)
        {
          double Pshock = density * (intercept_EOS*1000 + slope_EOS*particle_velocity*1000)*particle_velocity*1000* std::pow((iso_core_radius/distance_to_impact),-1.84+2.61*std::log(velocity_impactor));
          double beta = density * std::pow(intercept_EOS*1000,2.0)/(2*slope_EOS);
          double f = -( Pshock / beta ) *std::pow(( 1 - std::sqrt( 2*Pshock/beta + 1) ), -1);
          double temperature_impact = std::pow(specific_heat,-1) * ((Pshock/(2*density*slope_EOS))  *(1-std::pow(f,-1)) - (f-std::log(f)-1)*std::pow(intercept_EOS*1000/slope_EOS, 2.0)); 
        }
      
      return temperature_impact;

      // else if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::cartesian)
      //   {
      //     Point<dim> impact_point;
      //     if (dim==2)
      //       {
      //         impact_point(0) = impact_position_one,
      //         impact_point(1) = Plugins::get_plugin_as_type<const GeometryModel::Box<dim>>(this->get_geometry_model()).get_extents()[dim-1];
      //       }
      //     else if (dim==3)
      //       {
      //         impact_point(0) = impact_position_one,
      //         impact_point(1) = impact_position_two,
      //         impact_point(2) = Plugins::get_plugin_as_type<const GeometryModel::Box<dim>>(this->get_geometry_model()).get_extents()[dim-1];
      //       }
      //     else
      //       return 100;
          
      //     double distance_to_impact = std::abs(impact_point.distance(position));
      //           //third, look up material properties and calculate the impact temperature for each points
      //     typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, this->n_compositional_fields());
      //     typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, this->n_compositional_fields());
      //     if (distance_to_impact<iso_core_radius)
      //     {
      //       double Pshock = density * (intercept_EOS*1000 + slope_EOS*particle_velocity*1000)*particle_velocity*1000;
      //       double beta = density * std::pow(intercept_EOS*1000,2.0)/(2*slope_EOS);
      //       double f = -( Pshock / beta ) *std::pow(( 1 - std::sqrt( 2*Pshock/beta + 1) ), -1);
      //       double temperature_impact = std::pow(specific_heat,-1) * ((Pshock/(2*density*slope_EOS))  *(1-std::pow(f,-1)) - (f-std::log(f)-1)*std::pow(intercept_EOS*1000/slope_EOS, 2.0)); 
      //       return temperature_impact;
      //     }
      //     else if (distance_to_impact>=iso_core_radius)
      //     {
      //       double Pshock = density * (intercept_EOS*1000 + slope_EOS*particle_velocity*1000)*particle_velocity*1000* std::pow((iso_core_radius/distance_to_impact),-1.84+2.61*std::log(velocity_impactor));
      //       double beta = density * std::pow(intercept_EOS*1000,2.0)/(2*slope_EOS);
      //       double f = -( Pshock / beta ) *std::pow(( 1 - std::sqrt( 2*Pshock/beta + 1) ), -1);
      //       double temperature_impact = std::pow(specific_heat,-1) * ((Pshock/(2*out.densities[0]*slope_EOS))  *(1-std::pow(f,-1)) - (f-std::log(f)-1)*std::pow(intercept_EOS*1000/slope_EOS, 2.0)); 
      //       return temperature_impact;
      //     }
      //     else
      //       return 100;
      //   }
      //   else
      //     return 100;
      //calculate sum and (return) limited temperature
    }
  


 //declare and parse parameters
    template <int dim>
    void
    Impact<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Impact");
        {
          prm.declare_entry ("Impact Position One", "0.0",
                             Patterns::Double (0),
                             "The position coordinates of the impactor "
                             "This value represents x for cartesian, and phi for spherical."
                             "Value: 0~x_max or 0~2pi.");
          prm.declare_entry ("Impact Position Two", "0.0",
                             Patterns::Double (0),
                             "The position coordinates of the impactor "
                             "This value represents y for cartesian, and theta for spherical."
                             "Value: 0~y_max or 0~pi");
          prm.declare_entry ("Ri", "100.0",
                             Patterns::Double (0),
                             "Radius of the impactor." 
                             "Units: km");
          prm.declare_entry ("Vi", "15.0",
                             Patterns::Double (0),
                             "Velocity of the impactor."
                             "Units: km/s");
          prm.declare_entry ("Uc", "7.5",
                             Patterns::Double (0),
                             "The shock-induced particle velocity within the isobaric core"
                             "Under the assumption that the target and projectile materials are the same, Uc=(Vi)/2"
                             "Units: km/s");
          prm.declare_entry ("C", "5.2",
                             Patterns::Double (0),
                             "Intercept of the linear Hugoniot shock equation of state (EOS) "
                             "C is also the acoustic velocity in the target material "
                             "Units: km/s");
          prm.declare_entry ("S", "1.5",
                             Patterns::Double (0),
                             "Slope of the linear Hugoniot shock equation of state (EOS) ");
          prm.declare_entry ("Reference Density", "3400",
                             Patterns::Double (0),
                             "This value is only a reference density for calculating impact heating. "
                             "Units: kg/m3");
          prm.declare_entry ("Reference Specific Heat", "1200",
                             Patterns::Double (0),
                             "This value is only for calculating impact heating. "
                             "Units: J/(kg*K)");
          // prm.declare_entry ("A1", "1085.7",
          //                    Patterns::Double (),
          //                    "Constant parameter of the solidus."
          //                    "Units: \\si{\\degreeCelsius}.");
          // prm.declare_entry ("A2", "1.329e-7",
          //                    Patterns::Double (),
          //                    "Prefactor of the linear pressure term of the solidus."
          //                    "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          // prm.declare_entry ("A3", "-5.1e-18",
          //                    Patterns::Double (),
          //                    "Prefactor of the quadratic pressure term of the solidus."
          //                    "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
          // prm.declare_entry ("B1", "1475.0",
          //                    Patterns::Double (),
          //                    "Constant parameter of the liquidus."
          //                    "Units: \\si{\\degreeCelsius}.");
          // prm.declare_entry ("B2", "8.0e-8",
          //                    Patterns::Double (),
          //                    "Prefactor of the linear pressure term of the liquidus."
          //                    "Units: \\si{\\degreeCelsius\\per\\pascal}.");
          // prm.declare_entry ("B3", "-3.2e-18",
          //                    Patterns::Double (),
          //                    "Prefactor of the quadratic pressure term of the liquidus."
          //                    "Units: \\si{\\degreeCelsius\\per\\pascal\\squared}.");
        
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    Impact<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Impact");
        {
          impact_position_one = prm.get_double ("Impact Position One");
          impact_position_two = prm.get_double ("Impact Position Two");
          radius_impactor = prm.get_double ("Ri");
          velocity_impactor = prm.get_double ("Vi");
          particle_velocity = prm.get_double ("Uc");
          intercept_EOS = prm.get_double ("C");
          slope_EOS = prm.get_double ("S");
          reference_density_value = prm.get_double ("Reference Density");
          reference_specific_heat_value = prm.get_double ("Reference Specific Heat");
          // A1              = prm.get_double ("A1");
          // A2              = prm.get_double ("A2");
          // A3              = prm.get_double ("A3");
          // B1              = prm.get_double ("B1");
          // B2              = prm.get_double ("B2");
          // B3              = prm.get_double ("B3");
          // solidus\liquidus\none
          // function:initial temperature
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Impact,
                                              "impact",
                                              "Impact-induced shock heating")
  }
}
